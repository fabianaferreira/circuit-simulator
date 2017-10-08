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

#define versao "1.0"
#define MAX_LINHA 80
#define MAX_TIPO_FONTE  5
#define MAX_NOME 11
#define MAX_ELEM 50
#define MAX_NOS 50
#define TOLG 1e-9
//#define DEBUG
#define PI acos(-1.0)
#define NOME_ARQUIVO_SAIDA "saida_simulacao.tab"

typedef struct sine /* STRUCT SIN */
{
  double nivel_dc,
         amplitude,
         freq,
         atraso,
         amortecimento,
         defasagem;
  unsigned int ciclos;
} sine;

typedef struct dc /* STRUCT DC */
{
  double valor;
} dc;

typedef struct pulse /* STRUCT PULSE */
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
int fontes_variantes[MAX_ELEM];

int
  numeroElementos, /* Elementos */
  numeroVariaveis, /* Variaveis */
  numeroNos, /* Nos */
  i,j,k;

char
/* Foram colocados limites nos formatos de leitura para alguma protecao
   contra excesso de caracteres nestas variaveis */
  nomeArquivo[MAX_LINHA+1],
  tipo,
  na[MAX_NOME],nb[MAX_NOME],nc[MAX_NOME],nd[MAX_NOME],
  lista[MAX_NOS+1][MAX_NOME+2], /*Tem que caber jx antes do nome */
  txt[MAX_LINHA+1],
  *p;
FILE *arquivo;

double
  g,
  Yn[MAX_NOS+1][MAX_NOS+2],
  YnAnterior[MAX_NOS+1][MAX_NOS+2];

/*  Rotina para Resolucao de sistema de equacoes lineares.
   Metodo de Gauss-Jordan com condensacao pivotal */
int ResolverSistema(void)
{
  int i,j,l, a;
  double t, p;

  for (i=1; i<=numeroVariaveis; i++) {
    t=0.0;
    a=i;
    for (l=i; l<=numeroVariaveis; l++) {
      if (fabs(Yn[l][i])>fabs(t)) {
	a=l;
	t=Yn[l][i];
      }
    }
    if (i!=a) {
      for (l=1; l<=numeroVariaveis+1; l++) {
	p=Yn[i][l];
	Yn[i][l]=Yn[a][l];
	Yn[a][l]=p;
      }
    }
    if (fabs(t)<TOLG) {
      printf("Sistema singular\n");
      return 1;
    }
    for (j=numeroVariaveis+1; j>0; j--) {  /* Basta j>i em vez de j>0 */
      Yn[i][j]/= t;
      p=Yn[i][j];
      if (p!=0)  /* Evita operacoes com zero */
        for (l=1; l<=numeroVariaveis; l++) {
	  if (l!=i)
	    Yn[l][j]-=Yn[l][i]*p;
        }
    }
  }
  return 0;
}

/* Rotina que conta os nos e atribui numeros a eles */
int NumerarNo(char *nome)
{
  int i,achou;

  i=0; achou=0;
  while (!achou && i<=numeroVariaveis)
    if (!(achou=!strcmp(nome,lista[i]))) i++;
  if (!achou) {
    if (numeroVariaveis==MAX_NOS) {
      printf("O programa so aceita ate %d nos\n",numeroVariaveis);
      exit(1);
    }
    numeroVariaveis++;
    strcpy(lista[numeroVariaveis],nome);
    return numeroVariaveis; /* novo no */
  }
  else {
    return i; /* no ja conhecido */
  }
}

void ArmazenarResultadoAnterior()
{
  for (i=0; i<=numeroVariaveis; i++)
  {
    for (j=0; j<=numeroVariaveis+1; j++)
      YnAnterior[i][j] = Yn[i][j];
  }
}

void ZerarSistema()
{
  for (i=0; i<=numeroVariaveis; i++)
  {
    for (j=0; j<=numeroVariaveis+1; j++)
      Yn[i][j]=0;
  }
}

/* Rotina que monta as estampas do circuito */
void MontarEstampas(double tempo, int pontoOp)
{
  unsigned ciclos, ciclos_passados;
  elemento fonte_atual;
  double amplitude,
         amplitude1,
         amplitude2,
         nivel_dc,
         atraso,
         freq,
         defasagem,
         amortecimento,
         tempo_subida,
         tempo_descida,
         tempo_ligada,
         periodo,
         tempo_normalizado;
    double coefAng,
           coefLin,
           t1,
           t2;
	for (i=1; i<=numeroElementos; i++)
	{
		tipo=netlist[i].nome[0];
		if (tipo=='R')
		{
			g = 1/netlist[i].valor;
			Yn[netlist[i].a][netlist[i].a] += g;
			Yn[netlist[i].b][netlist[i].b] += g;
			Yn[netlist[i].a][netlist[i].b] -= g;
			Yn[netlist[i].b][netlist[i].a] -= g;
		}
		else if (tipo=='L')
		{
			// Vira um R de 1nOhms
			g = 1/1e-9;
			Yn[netlist[i].a][netlist[i].x] += 1;
			Yn[netlist[i].b][netlist[i].x] -= 1;
			Yn[netlist[i].x][netlist[i].a] -= 1;
			Yn[netlist[i].x][netlist[i].b] += 1;
			Yn[netlist[i].x][netlist[i].x] += 1/g;
		}
		else if (tipo=='C')
		{
			// Vira um R de 1GOhms
			g = 1/1e9;
			Yn[netlist[i].a][netlist[i].a] += g;
			Yn[netlist[i].b][netlist[i].b] += g;
			Yn[netlist[i].a][netlist[i].b] -= g;
			Yn[netlist[i].b][netlist[i].a] -= g;
		}
		else if (tipo=='G')
		{
			g = netlist[i].valor;
			Yn[netlist[i].a][netlist[i].c] += g;
			Yn[netlist[i].b][netlist[i].d] += g;
			Yn[netlist[i].a][netlist[i].d] -= g;
			Yn[netlist[i].b][netlist[i].c] -= g;
		}
    else if (tipo=='I')
    {
      fonte_atual = netlist[i];
      if (strcmp(fonte_atual.tipo_fonte, "SIN") == 0)
      {
        amplitude = fonte_atual.fonte_seno.amplitude;
        freq = fonte_atual.fonte_seno.freq;
        atraso = fonte_atual.fonte_seno.atraso;
        defasagem = fonte_atual.fonte_seno.defasagem;
        nivel_dc = fonte_atual.fonte_seno.nivel_dc;
        amortecimento = fonte_atual.fonte_seno.amortecimento;
        ciclos = fonte_atual.fonte_seno.ciclos;
        fonte_atual.valor = nivel_dc +
                            amplitude*(exp(-amortecimento*(tempo - atraso)))*(sin(2*PI*freq*(tempo - atraso) + (PI*defasagem)/180));
      }
      else if (strcmp(netlist[i].tipo_fonte, "PULSE") == 0)
      {
        amplitude1 = fonte_atual.fonte_pulso.amplitude1;
        amplitude2 = fonte_atual.fonte_pulso.amplitude2;
        atraso = fonte_atual.fonte_pulso.atraso;
        tempo_subida = fonte_atual.fonte_pulso.tempo_subida;
        tempo_descida = fonte_atual.fonte_pulso.tempo_descida;
        ciclos = fonte_atual.fonte_pulso.ciclos;
        periodo = fonte_atual.fonte_pulso.periodo;
        tempo_ligada = fonte_atual.fonte_pulso.tempo_ligada;
        ciclos_passados = unsigned(tempo/periodo);

        /*Tratando descontinuidades*/
        if (tempo_subida == 0)
          tempo_descida = passo_simulacao;
        if (tempo_descida == 0)
          tempo_subida = passo_simulacao;

        tempo_normalizado = tempo - periodo*ciclos_passados;
        /*Falta o que fazer quando tá terminando*/
        /*Achando o valor da fonte no tempo atual*/
        if (ciclos_passados >= ciclos)
          fonte_atual.valor = amplitude1;
        else if (tempo_normalizado <= atraso)
          fonte_atual.valor = amplitude1;
        else if (tempo_normalizado <= tempo_subida + atraso)
        {
          /*Achando a equacao da reta de subida*/
          t1 = atraso;
          t2 = atraso + tempo_subida;
          coefAng = (amplitude2 - amplitude1)/(t2 - t1);
          coefLin = amplitude1 - coefAng*t1;
          fonte_atual.valor = coefAng*tempo + coefLin; /*????????????*/
        }
        else if (tempo_normalizado <= atraso + tempo_subida + tempo_ligada)
          fonte_atual.valor = amplitude2;
        else if (tempo_normalizado <= periodo)
        {
          /*Achando a equacao da reta de descida*/
          t1 = atraso + tempo_subida + tempo_ligada;
          t2 = t1 + tempo_descida;
          coefAng = (amplitude1 - amplitude2)/(t1 - t2);
          coefLin = amplitude1 - coefAng*t1;
          fonte_atual.valor = coefAng*tempo + coefLin;
        }
      }
      else
      {
        fonte_atual.valor = fonte_atual.fonte_dc.valor;
      }
      g=fonte_atual.valor;
      Yn[fonte_atual.a][numeroVariaveis+1]-=g;
      Yn[fonte_atual.b][numeroVariaveis+1]+=g;
    }

    /*Monta a estampa apenas se a fonte for DC*/
    else if (tipo=='V')
    {
      fonte_atual = netlist[i];
      if (strcmp(fonte_atual.tipo_fonte, "SIN") == 0)
      {
        amplitude = fonte_atual.fonte_seno.amplitude;
        freq = fonte_atual.fonte_seno.freq;
        atraso = fonte_atual.fonte_seno.atraso;
        defasagem = fonte_atual.fonte_seno.defasagem;
        nivel_dc = fonte_atual.fonte_seno.nivel_dc;
        amortecimento = fonte_atual.fonte_seno.amortecimento;
        ciclos = fonte_atual.fonte_seno.ciclos;
        fonte_atual.valor = nivel_dc +
                            amplitude*(exp(-amortecimento*(tempo - atraso)))*(sin(2*PI*freq*(tempo - atraso) + (PI*defasagem)/180));
      }
      else if (strcmp(fonte_atual.tipo_fonte, "PULSE") == 0)
      {
        amplitude1 = fonte_atual.fonte_pulso.amplitude1;
        amplitude2 = fonte_atual.fonte_pulso.amplitude2;
        atraso = fonte_atual.fonte_pulso.atraso;
        tempo_subida = fonte_atual.fonte_pulso.tempo_subida;
        tempo_descida = fonte_atual.fonte_pulso.tempo_descida;
        ciclos = fonte_atual.fonte_pulso.ciclos;
        periodo = fonte_atual.fonte_pulso.periodo;
        tempo_ligada = fonte_atual.fonte_pulso.tempo_ligada;
        ciclos_passados = unsigned(tempo/periodo);

        /*Tratando descontinuidades*/
        if (tempo_subida == 0)
          tempo_descida = passo_simulacao;
        if (tempo_descida == 0)
          tempo_subida = passo_simulacao;

        tempo_normalizado = tempo - periodo*ciclos_passados;
        /*Falta o que fazer quando tá terminando*/
        /*Achando o valor da fonte no tempo atual*/
        if (ciclos_passados >= ciclos)
          fonte_atual.valor = amplitude1;
        else if (tempo_normalizado <= atraso)
          fonte_atual.valor = amplitude1;
        else if (tempo_normalizado <= tempo_subida + atraso)
        {
          /*Achando a equacao da reta de subida*/
          t1 = atraso;
          t2 = atraso + tempo_subida;
          coefAng = (amplitude2 - amplitude1)/(t2 - t1);
          coefLin = amplitude1 - coefAng*t1;
          fonte_atual.valor = coefAng*tempo + coefLin; /*????????????*/
        }
        else if (tempo_normalizado <= atraso + tempo_subida + tempo_ligada)
          fonte_atual.valor = amplitude2;
        else if (tempo_normalizado <= periodo)
        {
          /*Achando a equacao da reta de descida*/
          t1 = atraso + tempo_subida + tempo_ligada;
          t2 = t1 + tempo_descida;
          coefAng = (amplitude1 - amplitude2)/(t1 - t2);
          coefLin = amplitude1 - coefAng*t1;
          fonte_atual.valor = coefAng*tempo + coefLin;
        }
      }
      else
      {
        fonte_atual.valor = fonte_atual.fonte_dc.valor;
      }
      g=fonte_atual.valor;
      Yn[fonte_atual.a][fonte_atual.x]+=1;
      Yn[fonte_atual.b][fonte_atual.x]-=1;
      Yn[fonte_atual.x][fonte_atual.a]-=1;
      Yn[fonte_atual.x][fonte_atual.b]+=1;
      Yn[fonte_atual.x][numeroVariaveis+1]-=g;
    }
		else if (tipo=='E')
		{
			g = netlist[i].valor;
			Yn[netlist[i].a][netlist[i].x] += 1;
			Yn[netlist[i].b][netlist[i].x] -= 1;
			Yn[netlist[i].x][netlist[i].a] -= 1;
			Yn[netlist[i].x][netlist[i].b] += 1;
			Yn[netlist[i].x][netlist[i].c] += g;
			Yn[netlist[i].x][netlist[i].d] -= g;
		}
		else if (tipo=='F')
		{
			g = netlist[i].valor;
			Yn[netlist[i].a][netlist[i].x] += g;
			Yn[netlist[i].b][netlist[i].x] -= g;
			Yn[netlist[i].c][netlist[i].x] += 1;
			Yn[netlist[i].d][netlist[i].x] -= 1;
			Yn[netlist[i].x][netlist[i].c] -= 1;
			Yn[netlist[i].x][netlist[i].d] += 1;
		}
		else if (tipo=='H')
		{
			g = netlist[i].valor;
			Yn[netlist[i].a][netlist[i].y] += 1;
			Yn[netlist[i].b][netlist[i].y] -= 1;
			Yn[netlist[i].c][netlist[i].x] += 1;
			Yn[netlist[i].d][netlist[i].x] -= 1;
			Yn[netlist[i].y][netlist[i].a] -= 1;
			Yn[netlist[i].y][netlist[i].b] += 1;
			Yn[netlist[i].x][netlist[i].c] -= 1;
			Yn[netlist[i].x][netlist[i].d] += 1;
			Yn[netlist[i].y][netlist[i].x] += g;
		}
		else if (tipo=='O')
		{
			Yn[netlist[i].a][netlist[i].x] += 1;
			Yn[netlist[i].b][netlist[i].x] -= 1;
			Yn[netlist[i].x][netlist[i].c] += 1;
			Yn[netlist[i].x][netlist[i].d] -= 1;
		}
#ifdef DEBUG
		/* Opcional: Mostra o sistema apos a montagem da estampa */
		printf("Sistema apos a estampa de %s\n",netlist[i].nome);
		for (k=1; k<=numeroVariaveis; k++)
		{
			for (j=1; j<=numeroVariaveis+1; j++)
				if (Yn[k][j]!=0)
					printf("%+4.3f ",Yn[k][j]);
				else printf(" ..... ");
			printf("\n");
		}
		getch();
#endif
	}
}

int main(void)
{
  unsigned contador_fontes_variantes = 0;
  system("cls");
  printf("Programa demonstrativo de analise nodal modificada no dominio do tempo\n\n");
  printf("Desenvolvido por: - Fabiana Ferreira Fonseca \n");
  printf("\t\t  - Leonardo Barreto Alves \n");
  printf("\t\t  - Vinicius dos Santos Mello \n\n");
  printf("\t\t  __/\\  /\\  /\\  /\\__\n");
  printf("\t\t      \\/  \\/  \\/    \n");
  printf("Versao %s\n",versao);
 denovo:
  /* Leitura do netlist */
  numeroElementos=0; numeroVariaveis=0; strcpy(lista[0],"0");
  printf("Nome do arquivo com o netlist (ex: mna.net): ");
  scanf("%50s",nomeArquivo);
  arquivo=fopen(nomeArquivo,"r");
  if (arquivo==0)
  {
    printf("Arquivo %s inexistente\n",nomeArquivo);
    goto denovo;
  }
  printf("Lendo netlist:\n");
  fgets(txt,MAX_LINHA,arquivo);
  printf("Titulo: %s",txt);
  while (fgets(txt,MAX_LINHA,arquivo))
  {
    numeroElementos++; /* Nao usa o netlist[0] */
    if (numeroElementos>MAX_ELEM)
    {
      printf("O programa so aceita ate %d elementos\n",MAX_ELEM);
      exit(1);
    }
    txt[0]=toupper(txt[0]);
    tipo=txt[0];
    sscanf(txt,"%10s",netlist[numeroElementos].nome);
    p=txt+strlen(netlist[numeroElementos].nome); /* Inicio dos parametros */
    /* O que e lido depende do tipo */

    /*RESISTOR*/
    if (tipo=='R')
    {
      sscanf(p,"%10s%10s%lg",na,nb,&netlist[numeroElementos].valor);
      printf("%s %s %s %g\n",netlist[numeroElementos].nome,na,nb,netlist[numeroElementos].valor);
      netlist[numeroElementos].a=NumerarNo(na);
      netlist[numeroElementos].b=NumerarNo(nb);
    }
    /*CAPACITOR*/
    else if (tipo=='C')
    {
      sscanf(p,"%10s%10s%lg",na,nb,&netlist[numeroElementos].valor);
      printf("%s %s %s %g\n",netlist[numeroElementos].nome,na,nb,netlist[numeroElementos].valor);
      netlist[numeroElementos].a=NumerarNo(na);
      netlist[numeroElementos].b=NumerarNo(nb);
    }
    /*INDUTOR*/
    else if (tipo=='L')
    {
      sscanf(p,"%10s%10s%lg",na,nb,&netlist[numeroElementos].valor);
      printf("%s %s %s %g\n",netlist[numeroElementos].nome,na,nb,netlist[numeroElementos].valor);
      netlist[numeroElementos].a=NumerarNo(na);
      netlist[numeroElementos].b=NumerarNo(nb);
    }

    /*
      No caso de fontes de tensao ou de corrente, precisamos identificar qual tipo de fonte
      para que possamo pegar os valores corretos, como amplitude e frequencia, no caso de fontes
      não constantes. Para tanto, sao usadas tres structs, que passam a ser atributos do struct elemento
    */

    else if (tipo == 'I' || tipo == 'V')
    {
      sscanf(p,"%10s%10s%5s",na,nb,netlist[numeroElementos].tipo_fonte);

      if (strcmp(netlist[numeroElementos].tipo_fonte, "DC") == 0)
      {
        /*valor*/
        sscanf(p, "%*10s%*10s%*5s%lg", &netlist[numeroElementos].fonte_dc.valor);
      }
      else if (strcmp(netlist[numeroElementos].tipo_fonte, "SIN") == 0)
      {
        /*nivel_dc, amplitude, freq, atraso, amortecimento, defasagem, ciclos;*/
        sscanf(p, "%*10s%*10s%*5s%lg%lg%lg%lg%lg%lg%i", &netlist[numeroElementos].fonte_seno.nivel_dc,
                &netlist[numeroElementos].fonte_seno.amplitude, &netlist[numeroElementos].fonte_seno.freq,
                &netlist[numeroElementos].fonte_seno.atraso, &netlist[numeroElementos].fonte_seno.amortecimento,
                &netlist[numeroElementos].fonte_seno.defasagem, &netlist[numeroElementos].fonte_seno.ciclos);
        fontes_variantes[contador_fontes_variantes] = numeroElementos;
        contador_fontes_variantes++;
      }
      else if (strcmp(netlist[numeroElementos].tipo_fonte, "PULSE") == 0)
      {
        /*amplitude1, amplitude2, atraso, tempo_subida,tempo_descida,tempo_ligada, periodo, ciclos;*/
        sscanf(p, "%*10s%*10s%*5s%lg%lg%lg%lg%lg%lg%lg%i", &netlist[numeroElementos].fonte_pulso.amplitude1,
                &netlist[numeroElementos].fonte_pulso.amplitude2, &netlist[numeroElementos].fonte_pulso.atraso,
                &netlist[numeroElementos].fonte_pulso.tempo_subida, &netlist[numeroElementos].fonte_pulso.tempo_descida,
                &netlist[numeroElementos].fonte_pulso.tempo_ligada, &netlist[numeroElementos].fonte_pulso.periodo,
                &netlist[numeroElementos].fonte_pulso.ciclos);
        fontes_variantes[contador_fontes_variantes] = numeroElementos;
        contador_fontes_variantes++;
      }
      netlist[numeroElementos].a=NumerarNo(na);
      netlist[numeroElementos].b=NumerarNo(nb);
    }

    /*FONTES CONTROLADAS*/
    else if (tipo=='G' || tipo=='E' || tipo=='F' || tipo=='H')
    {
      sscanf(p,"%10s%10s%10s%10s%lg",na,nb,nc,nd,&netlist[numeroElementos].valor);
      printf("%s %s %s %s %s %g\n",netlist[numeroElementos].nome,na,nb,nc,nd,netlist[numeroElementos].valor);
      netlist[numeroElementos].a=NumerarNo(na);
      netlist[numeroElementos].b=NumerarNo(nb);
      netlist[numeroElementos].c=NumerarNo(nc);
      netlist[numeroElementos].d=NumerarNo(nd);
    }

    /*AMPLIFICADOR OPERACIONAL IDEAL*/
    else if (tipo=='O')
    {
      sscanf(p,"%10s%10s%10s%10s",na,nb,nc,nd);
      printf("%s %s %s %s %s\n",netlist[numeroElementos].nome,na,nb,nc,nd);
      netlist[numeroElementos].a=NumerarNo(na);
      netlist[numeroElementos].b=NumerarNo(nb);
      netlist[numeroElementos].c=NumerarNo(nc);
      netlist[numeroElementos].d=NumerarNo(nd);
    }
    /*TRANSFORMADOR IDEAL*/
    else if(tipo=='K')
    {
      sscanf(p,"%10s%10s%10s%10s%lg",na,nb,nc,nd,&netlist[numeroElementos].valor);
      printf("%s %s %s %s %s\n",netlist[numeroElementos].nome,na,nb,nc,nd);
      netlist[numeroElementos].a=NumerarNo(na);
      netlist[numeroElementos].b=NumerarNo(nb);
      netlist[numeroElementos].c=NumerarNo(nc);
      netlist[numeroElementos].d=NumerarNo(nd);
    }
    /* RESISTOR LINEAR POR PARTES */
    else if(tipo=='N')
    {

    }
    /* CHAVE */
    else if(tipo=='$')
    {

    }
    /* COMENTÁRIO */
    else if (tipo=='*')
    { /* Comentario comeca com "*" */
      printf("Comentario: %s",txt);
      numeroElementos--;
    }

    /*Atribuindo os valores dos passos de integracao e de escrita no arquivo de saida,
      alem do tempo total de simulacao definido no netlist*/
    else if (tipo == '.')
    {
      if (strcmp (netlist[numeroElementos].nome, ".TRAN") == 0)
      {
        sscanf(p, "%lg%lg%*10s%lg", &tempo_simulacao, &passo_simulacao, &passo_saida);
        printf("%lg %lg %lg\n", tempo_simulacao, passo_simulacao, passo_saida);
      }
      numeroElementos--;
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
  numeroNos=numeroVariaveis;
  for (i=1; i<=numeroElementos; i++)
  {
    tipo=netlist[i].nome[0];
    if (tipo=='V' || tipo=='E' || tipo=='F' || tipo=='O')
    {
      numeroVariaveis++;
      if (numeroVariaveis>MAX_NOS)
      {
        printf("As correntes extra excederam o numero de variaveis permitido (%d)\n",MAX_NOS);
        exit(1);
      }
      strcpy(lista[numeroVariaveis],"j"); /* Tem espaco para mais dois caracteres */
      strcat(lista[numeroVariaveis],netlist[i].nome);
      netlist[i].x=numeroVariaveis;
    }
    else if (tipo=='H')
    {
      numeroVariaveis=numeroVariaveis+2;
      if (numeroVariaveis>MAX_NOS)
      {
        printf("As correntes extra excederam o numero de variaveis permitido (%d)\n",MAX_NOS);
        exit(1);
      }
      strcpy(lista[numeroVariaveis-1],"jx"); strcat(lista[numeroVariaveis-1],netlist[i].nome);
      netlist[i].x=numeroVariaveis-1;
      strcpy(lista[numeroVariaveis],"jy"); strcat(lista[numeroVariaveis],netlist[i].nome);
      netlist[i].y=numeroVariaveis;
    }
  }
  getch();
  /* Lista tudo */
  printf("Variaveis internas: \n");
  for (i=0; i<=numeroVariaveis; i++)
    printf("%d -> %s\n",i,lista[i]);
  getch();
  printf("Netlist interno final\n");
  for (i=1; i<=numeroElementos; i++)
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
  printf("O circuito tem %d nos, %d variaveis e %d elementos\n",numeroNos,numeroVariaveis,numeroElementos);
  getch();

  FILE *arquivoSaida = fopen(NOME_ARQUIVO_SAIDA, "w");
  /*Escreve o header do arquivo de saida*/
  strcpy(txt,"");
  fprintf(arquivoSaida, "%s", "t");
  for (i=1; i<=numeroVariaveis; i++)
  {
    //if (i==numeroNos+1) strcpy(txt,"j");
    fprintf(arquivoSaida, " %s%s", txt, lista[i]);
  }
  fprintf(arquivoSaida, "\n");

  /*Analise no tempo*/
  for (tempo_atual = 0; tempo_atual < tempo_simulacao; tempo_atual += passo_simulacao)
  {
    ZerarSistema();
    MontarEstampas(tempo_atual, 0);

    /*Newton-Raphson para tempo atual com o ResolverSistema() varias vezes ate convergir*/

    /* Resolve o sistema */
    if (ResolverSistema())
    {
      getch();
      exit(0);
    }

    /*Escreve no arquivo de saida*/
    fprintf(arquivoSaida,"%lg", tempo_atual);
    for (i=1; i<=numeroVariaveis; i++)
    {
      fprintf(arquivoSaida," %lg", Yn[i][numeroVariaveis+1]);
    }
    fprintf(arquivoSaida,"\n");
  }/*for analise no tempo*/

  /*Fecha arquivo*/
  fclose(arquivoSaida);

#ifdef DEBUG
  /* Opcional: Mostra o sistema resolvido */
  printf("Sistema resolvido:\n");
  for (i=1; i<=numeroVariaveis; i++)
  {
      for (j=1; j<=numeroVariaveis+1; j++)
        if (Yn[i][j]!=0) printf("%+3.1f ",Yn[i][j]);
        else printf(" ... ");
      printf("\n");
    }
  getch();
#endif
  /* Mostra solucao */
  printf("Solucao:\n");
  strcpy(txt,"Tensao");
  for (i=1; i<=numeroVariaveis; i++)
  {
    if (i==numeroNos+1) strcpy(txt,"Corrente");
    printf("%s %s: %g\n",txt,lista[i],Yn[i][numeroVariaveis+1]);
  }
  getch();
  return 0;
}
