/*
Elementos aceitos e linhas do netlist:

Resistor:                     R<nome> <no+> <no-> <resistencia>
VCCS:                         G<nome> <io+> <io-> <vi+> <vi-> <transcondutancia>
VCVS:                         E<nome> <vo+> <vo-> <vi+> <vi-> <ganho de tensao>
CCCS:                         F<nome> <io+> <io-> <ii+> <ii-> <ganho de corrente>
CCVS:                         H<nome> <vo+> <vo-> <ii+> <ii-> <transresistencia>
Fonte I:                      I<nome> <no+> <no-> <parametros>
Fonte V:                      V<nome> <no+> <no-> <parametros>
Amp. op.:                     O<nome> <vo1> <vo2> <vi1> <vi2>
Resistor linear por partes:   N<nome> <no+> <no-> <4 pontos vi, ji>
Capacitor:                    C<nome> <no1> <no2> <capacitancia>
Indutor:                      L<nome> <no1> <no2> <indutancia>
Transformador ideal:          K<nome> <noa> <nob> <noc> <nod> <n>
Chave:                        $<nome> <noa> <nob> <noContc> <noContd> <gon> <goff> <vref>
Comentario:                   *<comentario>


Parametros:
Fonte contínua: DC <valor>
Fonte senoidal: SIN <nível contínuo> <amplitude> <frequência em Hz> <atraso> <amortecimento> <defasagem em graus> <número de ciclos>
Fonte pulsada: PULSE <amplitude 1> <amplitude 2> <atraso> <tempo de subida> <tempo de descida> <tempo ligada> <período> <número de ciclos>

As fontes F e H tem o ramo de entrada em curto
O amplificador operacional ideal tem a saida suspensa
Os nos podem ser nomes
*/

#include <stdio.h>
#include <conio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#define MAX_LINHA 80
#define MAX_TIPO_FONTE  5
#define MAX_NOME 11
#define MAX_ELEM 50
#define MAX_NOS 50
#define TOLG 1e-9
#define DEBUG

typedef struct sine /* CLASSE SIN */
{
  double nivel_dc,
         amplitude,
         freq,
         atraso,
         amortecimento,
         defasagem;
  unsigned int ciclos;
} sine;

typedef struct dc /* CLASSE DC */
{
  double valor;
} dc;

typedef struct pulse /* CLASSE PULSE */
{
  double amplitude1,
         amplitude2,
         atraso,
         tempo_subida,
         tempo_descida,
         tempo_ligada,
         periodo;
  unsigned int ciclos;
} pulse;

/*Elemento possui atributos de fontes, pois há tipos diferentes, com parametros diferentes*/
typedef struct elemento /* Elemento do netlist */
{
  char nome[MAX_NOME];
  double valor;
  char tipo_fonte[MAX_TIPO_FONTE];
  int a,b,c,d,x,y;
  sine fonte_seno;
  dc fonte_dc;
  pulse fonte_pulso;
} elemento;

/*As seguintes variaveis vao definir os passos e o tempo de simulacao a ser usado
  Como o passo a ser escrito no arquivo de saida pode nao ser o mesmo do passo da
  integracao, vamos definir os dois separadamente*/
double tempo_simulacao, passo_simulacao, passo_saida;
double tempo_atual;

elemento netlist[MAX_ELEM]; /* Netlist */

int
  ne, /* Elementos */
  nv, /* Variaveis */
  nn, /* Nos */
  i,j,k;

char
/* Foram colocados limites nos formatos de leitura para alguma protecao
   contra excesso de caracteres nestas variaveis */
  nomearquivo[MAX_LINHA+1],
  tipo,
  na[MAX_NOME],nb[MAX_NOME],nc[MAX_NOME],nd[MAX_NOME],
  lista[MAX_NOS+1][MAX_NOME+2], /*Tem que caber jx antes do nome */
  txt[MAX_LINHA+1],
  *p;
FILE *arquivo;

double
  g,
  Yn[MAX_NOS+1][MAX_NOS+2];

/* Resolucao de sistema de equacoes lineares.
   Metodo de Gauss-Jordan com condensacao pivotal */
int resolversistema(void)
{
  int i,j,l, a;
  double t, p;

  for (i=1; i<=nv; i++) {
    t=0.0;
    a=i;
    for (l=i; l<=nv; l++) {
      if (fabs(Yn[l][i])>fabs(t)) {
	a=l;
	t=Yn[l][i];
      }
    }
    if (i!=a) {
      for (l=1; l<=nv+1; l++) {
	p=Yn[i][l];
	Yn[i][l]=Yn[a][l];
	Yn[a][l]=p;
      }
    }
    if (fabs(t)<TOLG) {
      printf("Sistema singular\n");
      return 1;
    }
    for (j=nv+1; j>0; j--) {  /* Basta j>i em vez de j>0 */
      Yn[i][j]/= t;
      p=Yn[i][j];
      if (p!=0)  /* Evita operacoes com zero */
        for (l=1; l<=nv; l++) {
	  if (l!=i)
	    Yn[l][j]-=Yn[l][i]*p;
        }
    }
  }
  return 0;
}

/* Rotina que conta os nos e atribui numeros a eles */
int numero(char *nome)
{
  int i,achou;

  i=0; achou=0;
  while (!achou && i<=nv)
    if (!(achou=!strcmp(nome,lista[i]))) i++;
  if (!achou) {
    if (nv==MAX_NOS) {
      printf("O programa so aceita ate %d nos\n",nv);
      exit(1);
    }
    nv++;
    strcpy(lista[nv],nome);
    return nv; /* novo no */
  }
  else {
    return i; /* no ja conhecido */
  }
}

int main(void)
{
  //clrscr();
  printf("Programa demonstrativo de analise nodal modificada\n");
  printf("Por Antonio Carlos M. de Queiroz - acmq@coe.ufrj.br\n");
  printf("Versao %s\n",versao);
 denovo:
  /* Leitura do netlist */
  ne=0; nv=0; strcpy(lista[0],"0");
  printf("Nome do arquivo com o netlist (ex: mna.net): ");
  scanf("%50s",nomearquivo);
  arquivo=fopen(nomearquivo,"r");
  if (arquivo==0)
  {
    printf("Arquivo %s inexistente\n",nomearquivo);
    goto denovo;
  }
  printf("Lendo netlist:\n");
  fgets(txt,MAX_LINHA,arquivo);
  printf("Titulo: %s",txt);
  while (fgets(txt,MAX_LINHA,arquivo))
  {
    ne++; /* Nao usa o netlist[0] */
    if (ne>MAX_ELEM)
    {
      printf("O programa so aceita ate %d elementos\n",MAX_ELEM);
      exit(1);
    }
    txt[0]=toupper(txt[0]);
    tipo=txt[0];
    sscanf(txt,"%10s",netlist[ne].nome);
    p=txt+strlen(netlist[ne].nome); /* Inicio dos parametros */
    /* O que e lido depende do tipo */

    /*RESISTOR*/
    if (tipo=='R')
    {
      sscanf(p,"%10s%10s%lg",na,nb,&netlist[ne].valor);
      printf("%s %s %s %g\n",netlist[ne].nome,na,nb,netlist[ne].valor);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
    }

    /*
      No caso de fontes de tensao ou de corrente, precisamos identificar qual tipo de fonte
      para que possamo pegar os valores corretos, como amplitude e frequencia, no caso de fontes
      não constantes. Para tanto, sao usadas tres structs, que passam a ser atributos do struct elemento
    */
    if (tipo == 'I' || tipo == 'V')
    {
      sscanf(p,"%10s%10s%5s",na,nb,netlist[ne].tipo_fonte);

      if (strcmp(netlist[ne].tipo_fonte, "DC") == 0)
      {
        /*valor*/
        sscanf(p, "%*10s%*10s%*5s%lg", &netlist[ne].fonte_dc.valor);
      }

      else if (strcmp(netlist[ne].tipo_fonte, "SIN") == 0)
      {
        /*nivel_dc, amplitude, freq, atraso, amortecimento, defasagem, ciclos;*/
        sscanf(p, "%*10s%*10s%*5s%lg%lg%lg%lg%lg%lg%i", &netlist[ne].fonte_seno.nivel_dc,
                &netlist[ne].fonte_seno.amplitude, &netlist[ne].fonte_seno.freq,
                &netlist[ne].fonte_seno.atraso, &netlist[ne].fonte_seno.amortecimento,
                &netlist[ne].fonte_seno.defasagem, &netlist[ne].fonte_seno.ciclos);
      }

      else if (strcmp(netlist[ne].tipo_fonte, "PULSE") == 0)
      {
        /*amplitude1, amplitude2, atraso, tempo_subida,tempo_descida,tempo_ligada, periodo, ciclos;*/
        sscanf(p, "%*10s%*10s%*5s%lg%lg%lg%lg%lg%lg%lg%i", &netlist[ne].fonte_pulso.amplitude1,
                &netlist[ne].fonte_pulso.amplitude2, &netlist[ne].fonte_pulso.atraso,
                &netlist[ne].fonte_pulso.tempo_subida, &netlist[ne].fonte_pulso.tempo_descida,
                &netlist[ne].fonte_pulso.tempo_ligada, &netlist[ne].fonte_pulso.periodo,
                &netlist[ne].fonte_pulso.ciclos);
      }
    }

    /*FONTES CONTROLADAS*/
    else if (tipo=='G' || tipo=='E' || tipo=='F' || tipo=='H')
    {
      sscanf(p,"%10s%10s%10s%10s%lg",na,nb,nc,nd,&netlist[ne].valor);
      printf("%s %s %s %s %s %g\n",netlist[ne].nome,na,nb,nc,nd,netlist[ne].valor);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
      netlist[ne].c=numero(nc);
      netlist[ne].d=numero(nd);
    }

    /*AMPLIFICADOR OPERACIONAL IDEAL*/
    else if (tipo=='O')
    {
      sscanf(p,"%10s%10s%10s%10s",na,nb,nc,nd);
      printf("%s %s %s %s %s\n",netlist[ne].nome,na,nb,nc,nd);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
      netlist[ne].c=numero(nc);
      netlist[ne].d=numero(nd);
    }
    else if (tipo=='*')
    { /* Comentario comeca com "*" */
      printf("Comentario: %s",txt);
      ne--;
    }

    /*Atribuindo os valores dos passos de integracao e de escrita no arquivo de saida,
      alem do tempo total de simulacao definido no netlist*/
    else if (tipo == ".")
    {
      if (strcmp (netlist[ne].nome, ".TRAN") == 0)
      {
        sscanf(p, "%lg%lg%*10s%lg", &tempo_simulacao, &passo_simulacao, &passo_saida);
        printf("%lg %lg %lg\n", tempo_simulacao, passo_simulacao, passo_saida);
      }
      ne--;
    }

    else
    {
      printf("Elemento desconhecido: %s\n",txt);
      getch();
      exit(1);
    }
  }
  fclose(arquivo);

  /* Acrescenta variaveis de corrente acima dos nos, anotando no netlist */
  nn=nv;
  for (i=1; i<=ne; i++)
  {
    tipo=netlist[i].nome[0];
    if (tipo=='V' || tipo=='E' || tipo=='F' || tipo=='O')
    {
      nv++;
      if (nv>MAX_NOS)
      {
        printf("As correntes extra excederam o numero de variaveis permitido (%d)\n",MAX_NOS);
        exit(1);
      }
      strcpy(lista[nv],"j"); /* Tem espaco para mais dois caracteres */
      strcat(lista[nv],netlist[i].nome);
      netlist[i].x=nv;
    }
    else if (tipo=='H')
    {
      nv=nv+2;
      if (nv>MAX_NOS)
      {
        printf("As correntes extra excederam o numero de variaveis permitido (%d)\n",MAX_NOS);
        exit(1);
      }
      strcpy(lista[nv-1],"jx"); strcat(lista[nv-1],netlist[i].nome);
      netlist[i].x=nv-1;
      strcpy(lista[nv],"jy"); strcat(lista[nv],netlist[i].nome);
      netlist[i].y=nv;
    }
  }
  getch();
  /* Lista tudo */
  printf("Variaveis internas: \n");
  for (i=0; i<=nv; i++)
    printf("%d -> %s\n",i,lista[i]);
  getch();
  printf("Netlist interno final\n");
  for (i=1; i<=ne; i++)
  {
    tipo=netlist[i].nome[0];
    if (tipo=='R' || tipo=='I' || tipo=='V')
    {
      printf("%s %d %d %g\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].valor);
    }
    else if (tipo=='G' || tipo=='E' || tipo=='F' || tipo=='H')
    {
      printf("%s %d %d %d %d %g\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].c,netlist[i].d,netlist[i].valor);
    }
    else if (tipo=='O')
    {
      printf("%s %d %d %d %d\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].c,netlist[i].d);
    }
    if (tipo=='V' || tipo=='E' || tipo=='F' || tipo=='O')
      printf("Corrente jx: %d\n",netlist[i].x);
    else if (tipo=='H')
      printf("Correntes jx e jy: %d, %d\n",netlist[i].x,netlist[i].y);
  }
  getch();
  /* Monta o sistema nodal modificado */
  printf("O circuito tem %d nos, %d variaveis e %d elementos\n",nn,nv,ne);
  getch();
  /* Zera sistema */
  for (i=0; i<=nv; i++)
  {
    for (j=0; j<=nv+1; j++)
      Yn[i][j]=0;
  }
  /* Monta estampas */
  for (i=1; i<=ne; i++)
  {
    tipo=netlist[i].nome[0];
    if (tipo=='R')
    {
      g=1/netlist[i].valor;
      Yn[netlist[i].a][netlist[i].a]+=g;
      Yn[netlist[i].b][netlist[i].b]+=g;
      Yn[netlist[i].a][netlist[i].b]-=g;
      Yn[netlist[i].b][netlist[i].a]-=g;
    }
    else if (tipo=='G')
    {
      g=netlist[i].valor;
      Yn[netlist[i].a][netlist[i].c]+=g;
      Yn[netlist[i].b][netlist[i].d]+=g;
      Yn[netlist[i].a][netlist[i].d]-=g;
      Yn[netlist[i].b][netlist[i].c]-=g;
    }
    /*Monta a estampa apenas se a fonte for DC*/
    else if (tipo=='I' && (strcmp(netlist[i].tipo_fonte, "DC") == 0))
    {
      g=netlist[i].fonte_dc.valor;
      Yn[netlist[i].a][nv+1]-=g;
      Yn[netlist[i].b][nv+1]+=g;
    }

    /*Monta a estampa apenas se a fonte for DC*/
    else if (tipo=='V' && (strcmp(netlist[i].tipo_fonte, "DC") == 0))
    {
      Yn[netlist[i].a][netlist[i].x]+=1;
      Yn[netlist[i].b][netlist[i].x]-=1;
      Yn[netlist[i].x][netlist[i].a]-=1;
      Yn[netlist[i].x][netlist[i].b]+=1;
      Yn[netlist[i].x][nv+1]-=netlist[i].fonte_dc.valor;
    }
    else if (tipo=='E')
    {
      g=netlist[i].valor;
      Yn[netlist[i].a][netlist[i].x]+=1;
      Yn[netlist[i].b][netlist[i].x]-=1;
      Yn[netlist[i].x][netlist[i].a]-=1;
      Yn[netlist[i].x][netlist[i].b]+=1;
      Yn[netlist[i].x][netlist[i].c]+=g;
      Yn[netlist[i].x][netlist[i].d]-=g;
    }
    else if (tipo=='F')
    {
      g=netlist[i].valor;
      Yn[netlist[i].a][netlist[i].x]+=g;
      Yn[netlist[i].b][netlist[i].x]-=g;
      Yn[netlist[i].c][netlist[i].x]+=1;
      Yn[netlist[i].d][netlist[i].x]-=1;
      Yn[netlist[i].x][netlist[i].c]-=1;
      Yn[netlist[i].x][netlist[i].d]+=1;
    }
    else if (tipo=='H')
    {
      g=netlist[i].valor;
      Yn[netlist[i].a][netlist[i].y]+=1;
      Yn[netlist[i].b][netlist[i].y]-=1;
      Yn[netlist[i].c][netlist[i].x]+=1;
      Yn[netlist[i].d][netlist[i].x]-=1;
      Yn[netlist[i].y][netlist[i].a]-=1;
      Yn[netlist[i].y][netlist[i].b]+=1;
      Yn[netlist[i].x][netlist[i].c]-=1;
      Yn[netlist[i].x][netlist[i].d]+=1;
      Yn[netlist[i].y][netlist[i].x]+=g;
    }
    else if (tipo=='O')
    {
      Yn[netlist[i].a][netlist[i].x]+=1;
      Yn[netlist[i].b][netlist[i].x]-=1;
      Yn[netlist[i].x][netlist[i].c]+=1;
      Yn[netlist[i].x][netlist[i].d]-=1;
    }
    /*Agora, precisamos fazer a analise no tempo. Porem, as fontes SINE e PULSE tem valores
      que dependem do tempo de simulacao e do passo de integracao. Assim, as estampas das mesmas
      tem que ser montadas dentro de um loop que ira resolver o sistema para cada tempo diferente.
      Dessa forma, todas as estampas que nao dependem do tempo podem ser montadas antes, como esta
      sendo feita, mas, a partir daqui, iremos montar as estampas dependentes e resolver o sistema*/
    for (tempo_atual = 0; tempo_atual < tempo_simulacao; tempo_atual += passo_simulacao)
    {
      /*Vou definir uma funcao que vai alterar as estampas dos elementos variantes no tempo.
        Vai ser criado um vetor que, ao ler da netlist, vai conter os elementos que variam no
        tempo (no caso inicial, as fontes), para que, ao montar as estampas de tempos em tempos,
        o processo seja otimizado e não monte as estampas dos elementos que não sofrem variacao*/
    }

    #ifdef DEBUG
        /* Opcional: Mostra o sistema apos a montagem da estampa */
        printf("Sistema apos a estampa de %s\n",netlist[i].nome);
        for (k=1; k<=nv; k++)
        {
          for (j=1; j<=nv+1; j++)
            if (Yn[k][j]!=0) printf("%+3.1f ",Yn[k][j]);
            else printf(" ... ");
          printf("\n");
        }
        getch();
    #endif
  } /*end for monta estampas*/

  /* Resolve o sistema */
  if (resolversistema())
  {
    getch();
    exit;
  }
#ifdef DEBUG
  /* Opcional: Mostra o sistema resolvido */
  printf("Sistema resolvido:\n");
  for (i=1; i<=nv; i++)
  {
      for (j=1; j<=nv+1; j++)
        if (Yn[i][j]!=0) printf("%+3.1f ",Yn[i][j]);
        else printf(" ... ");
      printf("\n");
    }
  getch();
#endif
  /* Mostra solucao */
  printf("Solucao:\n");
  strcpy(txt,"Tensao");
  for (i=1; i<=nv; i++)
  {
    if (i==nn+1) strcpy(txt,"Corrente");
    printf("%s %s: %g\n",txt,lista[i],Yn[i][nv+1]);
  }
  getch();
  return 0;
}
