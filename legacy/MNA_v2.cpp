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
#include <time.h>

#define versao              "1.0"
#define MAX_LINHA           80
#define MAX_TIPO_FONTE      5
#define MAX_NOME            11
#define MAX_ELEM            50
#define MAX_NOS             50
#define TOLG                1e-9
#define PI acos(-1.0)
#define PO_CAPACITOR        1e9
#define PO_INDUTOR          1e-9
#define MAX_ERRO_NR         1e-6
//#define MAX_ERRO_GMIN       1e-5
#define X_ERRO              1
#define MAX_ITERACOES       100
#define MAX_INICIALIZACOES  10
#define GMIN_INICIAL        1.1
#define GMIN_MINIMA         1e-12
#define GMIN_MENOR_FATOR    1.1
#define NOME_ARQUIVO_SAIDA "saida_simulacao.tab"

//#define DEBUG
//#define DEBUG_GMIN

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

typedef struct chave
{
  double gon,
         goff,
         vLim;
} chave;

typedef struct resistor
{
  double v1, j1,
         v2, j2,
         v3, j3,
         v4, j4;
} resistor;

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
  chave chaveResistiva;
  resistor resistorPartes;
  double numeroEspiras;
  double jt0, vt0;
} elemento;

/*As seguintes variaveis vao definir os passos e o tempo de simulacao a ser usado
  Como o passo a ser escrito no arquivo de saida pode nao ser o mesmo do passo da
  integracao, vamos definir os dois separadamente*/
double tempo_simulacao, passo_simulacao;
//double tempo_atual;
double gminAtual = GMIN_INICIAL;
double fator = 10;
double g_anterior[MAX_ELEM];
double z_anterior [MAX_ELEM];

elemento netlist[MAX_ELEM]; /* Netlist */
int elementosVariantes[MAX_ELEM];
int elementosNaoLineares[MAX_ELEM];

int
  numeroElementos, /* Elementos */
  numeroVariaveis, /* Variaveis */
  numeroNos, /* Nos */
  i,j,k;

unsigned passosPorPonto = 1;
unsigned contadorElementosVariantes = 0;
unsigned contadorElementosNaoLineares = 0;
unsigned short temCapacitorOuIndutor = 0;

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
  YnNewtonRaphson[MAX_NOS+1][MAX_NOS+2],
  YnInvariantes[MAX_NOS+1][MAX_NOS+2],
  solucaoAnterior[MAX_NOS+2],
  newtonRaphsonAnterior[MAX_NOS+2],
  newtonRaphsonAnteriorConvergiu[MAX_NOS+2],
  erros[MAX_NOS+1];

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
  unsigned i;
  for (i=0; i<=numeroVariaveis+1; i++)
    solucaoAnterior[i] = Yn[i][numeroVariaveis+1]; /*pega a ultima coluna de Yn: solucao do sistema*/
}

void ArmazenarNRAnterior()
{
  unsigned i;
  for (i=0; i<=numeroVariaveis+1; i++)
    newtonRaphsonAnterior[i] = Yn[i][numeroVariaveis+1]; /*pega a ultima coluna de Yn: solucao do sistema*/
}

void ArmazenarUltimoNRConvergiu()
{
  unsigned i;
  for (i=0; i<=numeroVariaveis+1; i++)
    newtonRaphsonAnteriorConvergiu[i] = Yn[i][numeroVariaveis+1]; /*pega a ultima coluna de Yn: solucao do sistema*/
}

void ArmazenarResultadoNaoConvergencia()
{
  unsigned i;
  for (i=0; i<=numeroVariaveis+1; i++)
    solucaoAnterior[i] = newtonRaphsonAnteriorConvergiu[i]; /*pega a ultima coluna de Yn: solucao do sistema*/
}

void ZerarNRAnterior ()
{
  unsigned i;
  for (i=0; i<=numeroVariaveis+1; i++)
    newtonRaphsonAnterior[i] = 0;
}

void ArmazenarEstampasNewtonRaphson() /*vai armazenar a estampa do NR para adicionar o Gmin*/
{
  unsigned i,j;
  for (i=0; i<=numeroVariaveis; i++)
  {
    for (j=0; j<numeroVariaveis+1; j++)
      YnNewtonRaphson[i][j] = Yn[i][j];
  }
}

void CopiarEstampasNewtonRaphson ()
{
  unsigned i,j;
  for (i=0; i<=numeroVariaveis; i++)
  {
    for (j=0; j<numeroVariaveis+1; j++)
      Yn[i][j] = YnNewtonRaphson[i][j];
  }
}

void CopiarEstampaInvariante()
{
  unsigned i,j;
  for (i=0; i<=numeroVariaveis; i++)
  {
    for (j=0; j<=numeroVariaveis+1; j++)
      Yn[i][j] = YnInvariantes[i][j];
  }
}

void ZerarSistema()
{
  unsigned i,j;
  for (i=0; i<=numeroVariaveis; i++)
  {
    for (j=0; j<=numeroVariaveis+1; j++)
      Yn[i][j]=0;
  }
}

void MostrarSistema ()
{
  unsigned j,k;
  for (k=1; k<=numeroVariaveis; k++)
  {
  	for (j=1; j<=numeroVariaveis+1; j++)
  		if (Yn[k][j]!=0)
  			printf("%+4.3f ",Yn[k][j]);
  		else printf(" ..... ");
  	printf("\n");
  }
  getch();
}

void MostrarSolucaoAtual ()
{
  unsigned j;
  	for (j=1; j<=numeroVariaveis; j++)
  	 printf("%+4.3f ",Yn[j][numeroVariaveis+1]);
  	printf("\n");
  getch();
}

void LerNetlist (FILE *arquivo)
{
  contadorElementosVariantes = 0;
  contadorElementosNaoLineares = 0;
  printf("Lendo netlist:\n");
  fgets(txt,MAX_LINHA,arquivo);
  printf("Titulo: %s",txt);
  while (fgets(txt,MAX_LINHA,arquivo))
  {
    numeroElementos++; /* Nao usa o netlist[0] */
    if (numeroElementos > MAX_ELEM)
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
      elementosVariantes[contadorElementosVariantes] = numeroElementos;
      contadorElementosVariantes++;
      temCapacitorOuIndutor = 1;
    }
    /*INDUTOR*/
    else if (tipo=='L')
    {
      sscanf(p,"%10s%10s%lg",na,nb,&netlist[numeroElementos].valor);
      printf("%s %s %s %g\n",netlist[numeroElementos].nome,na,nb,netlist[numeroElementos].valor);
      netlist[numeroElementos].a=NumerarNo(na);
      netlist[numeroElementos].b=NumerarNo(nb);
      elementosVariantes[contadorElementosVariantes] = numeroElementos;
      contadorElementosVariantes++;
      temCapacitorOuIndutor = 1;
    }

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
        elementosVariantes[contadorElementosVariantes] = numeroElementos;
        contadorElementosVariantes++;
      }
      else if (strcmp(netlist[numeroElementos].tipo_fonte, "PULSE") == 0)
      {
        /*amplitude1, amplitude2, atraso, tempo_subida,tempo_descida,tempo_ligada, periodo, ciclos;*/
        sscanf(p, "%*10s%*10s%*5s%lg%lg%lg%lg%lg%lg%lg%i", &netlist[numeroElementos].fonte_pulso.amplitude1,
                &netlist[numeroElementos].fonte_pulso.amplitude2, &netlist[numeroElementos].fonte_pulso.atraso,
                &netlist[numeroElementos].fonte_pulso.tempo_subida, &netlist[numeroElementos].fonte_pulso.tempo_descida,
                &netlist[numeroElementos].fonte_pulso.tempo_ligada, &netlist[numeroElementos].fonte_pulso.periodo,
                &netlist[numeroElementos].fonte_pulso.ciclos);
        elementosVariantes[contadorElementosVariantes] = numeroElementos;
        contadorElementosVariantes++;
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
      //printf("Entrei no if do transformador\n");
      sscanf(p,"%10s%10s%10s%10s%lg",na,nb,nc,nd,&netlist[numeroElementos].numeroEspiras);
      printf("%s %s %s %s %s %lg\n",netlist[numeroElementos].nome,na,nb,nc,nd, netlist[numeroElementos].numeroEspiras);
      netlist[numeroElementos].a=NumerarNo(na);
      netlist[numeroElementos].b=NumerarNo(nb);
      netlist[numeroElementos].c=NumerarNo(nc);
      netlist[numeroElementos].d=NumerarNo(nd);
    }
    /* RESISTOR LINEAR POR PARTES */
    else if(tipo=='N')
    {
      sscanf(p,"%10s%10s%lg%lg%lg%lg%lg%lg%lg%lg",na,nb,&netlist[numeroElementos].resistorPartes.v1, &netlist[numeroElementos].resistorPartes.j1,
                                                        &netlist[numeroElementos].resistorPartes.v2, &netlist[numeroElementos].resistorPartes.j2,
                                                        &netlist[numeroElementos].resistorPartes.v3, &netlist[numeroElementos].resistorPartes.j3,
                                                        &netlist[numeroElementos].resistorPartes.v4, &netlist[numeroElementos].resistorPartes.j4);
      netlist[numeroElementos].a=NumerarNo(na);
      netlist[numeroElementos].b=NumerarNo(nb);
      elementosNaoLineares[contadorElementosNaoLineares] = numeroElementos;
      contadorElementosNaoLineares ++;
    }
    /* CHAVE */
    else if(tipo=='$')
    {
      printf("Achei as chaves\n");
      sscanf(p,"%10s%10s%10s%10s%lg%lg%lg",na,nb,nc,nd, &netlist[numeroElementos].chaveResistiva.gon, &netlist[numeroElementos].chaveResistiva.goff, &netlist[numeroElementos].chaveResistiva.vLim);
      printf("%s %s %s cntrl1:%s cntrl2:%s \n",netlist[numeroElementos].nome,na,nb,nc,nd);
      netlist[numeroElementos].a=NumerarNo(na);
      netlist[numeroElementos].b=NumerarNo(nb);
      netlist[numeroElementos].c=NumerarNo(nc);
      netlist[numeroElementos].d=NumerarNo(nd);

      elementosNaoLineares[contadorElementosNaoLineares] = numeroElementos;
      contadorElementosNaoLineares ++;
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
        sscanf(p, "%lg%lg%*10s%u", &tempo_simulacao, &passo_simulacao, &passosPorPonto);
        printf("%lg %lg %u\n", tempo_simulacao, passo_simulacao, passosPorPonto);
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
}/*LerNetlist*/

/*Pega apenas as fontes variantes, capacitores e indutores do netlist e monta a estampa de acordo com o tempo atual da simulacao*/
void MontarEstampasVariantes (double tempo, double passo_simulacao, unsigned pontoOperacao)
{
  //printf("Funcao MontarEstampasVariantes\n");
  unsigned contador, ciclos, ciclos_passados;
  elemento elementoVariante;
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
         tensaoAtual;

    if (pontoOperacao == 1)
      tempo = 0.0; /*fontes variantes calculadas em 0.0 para ponto de operacao*/

    CopiarEstampaInvariante();

    for (contador = 0; contador < contadorElementosVariantes; contador++)
    {
      elementoVariante = netlist[elementosVariantes[contador]];
      /*FONTE DE CORRENTE*/
      if (strcmp(elementoVariante.tipo_fonte, "SIN") == 0)
      {
        amplitude = elementoVariante.fonte_seno.amplitude;
        freq = elementoVariante.fonte_seno.freq;
        atraso = elementoVariante.fonte_seno.atraso;
        defasagem = elementoVariante.fonte_seno.defasagem;
        nivel_dc = elementoVariante.fonte_seno.nivel_dc;
        amortecimento = elementoVariante.fonte_seno.amortecimento;
        ciclos = elementoVariante.fonte_seno.ciclos;
        // ciclos_passados = freq*(tempo-atraso);
        //printf("Ciclos %u Ciclos passados %u\n", ciclos, ciclos_passados);

        if (tempo > atraso + ciclos/freq)
          tempo = atraso + ciclos/freq;

        if (tempo < atraso)
          elementoVariante.valor = nivel_dc + amplitude*(sin(PI*defasagem/180));
        else
          elementoVariante.valor = nivel_dc +
                                   amplitude*(exp(-amortecimento*(tempo - atraso)))*(sin(2*PI*freq*(tempo - atraso) + (PI*defasagem)/180));

        //printf("Valor fonte seno: %lg em t = %lg\n", elementoVariante.valor, tempo);
      }
      /*Tá dando merda, tem que tudo, inclusive os extremos*/
      else if (strcmp(elementoVariante.tipo_fonte, "PULSE") == 0)
      {
        //printf("Achou uma fonte pulso\n");
        amplitude1 = elementoVariante.fonte_pulso.amplitude1;
        amplitude2 = elementoVariante.fonte_pulso.amplitude2;
        atraso = elementoVariante.fonte_pulso.atraso;
        tempo_subida = elementoVariante.fonte_pulso.tempo_subida;
        tempo_descida = elementoVariante.fonte_pulso.tempo_descida;
        ciclos = elementoVariante.fonte_pulso.ciclos;
        periodo = elementoVariante.fonte_pulso.periodo;
        tempo_ligada = elementoVariante.fonte_pulso.tempo_ligada;

        /*Tratando descontinuidades*/
        if (tempo_subida == 0)
          tempo_subida = passo_simulacao;
        else
          tempo_subida = tempo_subida;

        if (tempo_descida == 0)
          tempo_descida = passo_simulacao;
        else
          tempo_descida = tempo_descida;

        // printf("Tempo de subida antes dos ifs: %lg\n", tempo_subida);

        double tempoReal = tempo - atraso;
        //printf("Tempo atual: %lg\n", tempo);
        tempoReal = fmod(tempoReal, periodo);
        // double tempoDesligada = periodo - (tempo_subida + tempo_descida + tempo_ligada);

        if (tempo <= atraso)
          elementoVariante.valor = amplitude1;

        else if (tempo >= (atraso + ciclos*periodo))
          elementoVariante.valor = amplitude1;

        else
        {
          if (tempoReal <= tempo_subida)
            elementoVariante.valor = ((amplitude2 - amplitude1)/tempo_subida)*tempoReal + amplitude1;

          else if (tempoReal <= tempo_subida + tempo_ligada + passo_simulacao)
            elementoVariante.valor = amplitude2;

          else if (tempoReal <= tempo_subida + tempo_ligada + tempo_descida)
            elementoVariante.valor = ((amplitude2 - amplitude1)/tempo_descida)*(tempoReal - (tempo_subida + tempo_ligada));

          else
            elementoVariante.valor = amplitude1;
        }
        // printf("Tempo de subida depois dos ifs: %lg\n", tempo_subida);

      } /*if pulse*/

      g=elementoVariante.valor;
      if (elementoVariante.nome[0] == 'I')
      {
        Yn[elementoVariante.a][numeroVariaveis+1]-=g;
        Yn[elementoVariante.b][numeroVariaveis+1]+=g;
      }

      else if (elementoVariante.nome[0] == 'V')
      {
        Yn[elementoVariante.a][elementoVariante.x]+=1;
        Yn[elementoVariante.b][elementoVariante.x]-=1;
        Yn[elementoVariante.x][elementoVariante.a]-=1;
        Yn[elementoVariante.x][elementoVariante.b]+=1;
        Yn[elementoVariante.x][numeroVariaveis+1]-=g;
      }

      /*VAO DEPENDER DO PASSO E DO VALOR ANTERIOR, SEJA PONTOOP OU RESULTADO DO NEWTON-RAPHSON*/

      else if (elementoVariante.nome[0]=='C') /*capacitor: elemento variante no tempo*/
      {
        if (pontoOperacao == 1) /*se é análise de ponto de operação*/
        {
          // Vira um R de 1GOhms
          g = 1/PO_CAPACITOR;
          //printf("Condutancia capacitor: %+4.20f\n", g);
    			Yn[elementoVariante.a][elementoVariante.a] += g;
    			Yn[elementoVariante.b][elementoVariante.b] += g;
    			Yn[elementoVariante.a][elementoVariante.b] -= g;
    			Yn[elementoVariante.b][elementoVariante.a] -= g;
        }
        else
        {
          g = (2*elementoVariante.valor)/passo_simulacao;
          Yn[elementoVariante.a][elementoVariante.a] += g;
    			Yn[elementoVariante.b][elementoVariante.b] += g;
    			Yn[elementoVariante.a][elementoVariante.b] -= g;
    			Yn[elementoVariante.b][elementoVariante.a] -= g;
          Yn[elementoVariante.a][numeroVariaveis+1]+= g*(elementoVariante.vt0) + elementoVariante.jt0;
          Yn[elementoVariante.b][numeroVariaveis+1]-= g*(elementoVariante.vt0) + elementoVariante.jt0;
        }
      }
      else if (elementoVariante.nome[0]=='L') /*indutor: elemento variante no tempo*/
      {
        if (pontoOperacao == 1) /*se é análise de ponto de operação*/
        {
          // Vira um R de 1nOhms
    			g = 1/PO_INDUTOR;
    			Yn[elementoVariante.a][elementoVariante.a] += g;
    			Yn[elementoVariante.b][elementoVariante.a] -= g;
    			Yn[elementoVariante.a][elementoVariante.b] -= g;
    			Yn[elementoVariante.b][elementoVariante.b] += g;
        }
        else
        {
          /*depende do passo e do valor anterior*/
          g = passo_simulacao/(2*elementoVariante.valor);
          Yn[elementoVariante.a][elementoVariante.a] += g;
          Yn[elementoVariante.a][elementoVariante.b] -= g;
          Yn[elementoVariante.b][elementoVariante.a] -= g;
          Yn[elementoVariante.b][elementoVariante.b] += g;

          Yn[elementoVariante.a][numeroVariaveis+1]-= (g*(elementoVariante.vt0) + elementoVariante.jt0);
          Yn[elementoVariante.b][numeroVariaveis+1]+= (g*(elementoVariante.vt0) + elementoVariante.jt0);
        }
      }
    }/*for*/

    #ifdef  DEBUG
      printf("Sistema apos montagem das estampas variantes no tempo. t = %g\n", tempo);
      MostrarSistema();
    #endif
}/*MontarEstampasVariantes*/

void InicializacaoRandomica()
{
	double valor;
  unsigned i;

	srand ((unsigned)time(NULL));

	//printf ("Inicializacao Randomica\n");
	for (i=1; i<=numeroVariaveis;i++)
	{
		if (erros[i] > MAX_ERRO_NR)
		{
			// printf ("Erro maior\n");
			if (i <= numeroNos)
			{
				valor = rand() % 20001;
				valor -= 10000;
				valor /= 1000.0;
				// printf ("i = %f\n", valor);
			}
			else
			{
				valor = rand() % 2001;
				valor -= 1000;
				valor /= 1000.0;
				// printf ("i = %f\n", valor);
			}
			newtonRaphsonAnterior[i] = valor;
		}
		else
			newtonRaphsonAnterior[i] = Yn[i][numeroVariaveis+1];
	}
	printf ("\n");
}

void MontarNewtonRaphson (double tempo, double passo_simulacao, unsigned int pontoOperacao) /*muda os elementos nao-lineares de acordo com o resultado anterior*/
{
  unsigned contador;
  elemento elementoNaoLinear;
  double tensaoAtual, g, z;
  /*gera a estampa novamente para ser alterada pelos elementos não-lineares*/
  ZerarSistema();
  CopiarEstampaInvariante();
  MontarEstampasVariantes(tempo, passo_simulacao, pontoOperacao);

  for (contador = 0; contador < contadorElementosNaoLineares; contador++)
  {
    elementoNaoLinear = netlist[elementosNaoLineares[contador]];
    if (elementoNaoLinear.nome[0] == '$') /*chave: elemento nao-linear*/
    {
      if (elementoNaoLinear.c == 0)
        newtonRaphsonAnterior[elementoNaoLinear.c] = 0;
      else if (elementoNaoLinear.d == 0)
        newtonRaphsonAnterior[elementoNaoLinear.d] = 0;

      tensaoAtual = newtonRaphsonAnterior[elementoNaoLinear.c] - newtonRaphsonAnterior[elementoNaoLinear.d];
      if (tensaoAtual < elementoNaoLinear.chaveResistiva.vLim)
      {
        z = 0;
        g = elementoNaoLinear.chaveResistiva.goff;
        Yn[elementoNaoLinear.a][elementoNaoLinear.a]+=g;
        Yn[elementoNaoLinear.b][elementoNaoLinear.b]+=g;
        Yn[elementoNaoLinear.a][elementoNaoLinear.b]-=g;
        Yn[elementoNaoLinear.b][elementoNaoLinear.a]-=g;
        g_anterior[contador]  = g;
        z_anterior[contador] = z;
      }
      else if (tensaoAtual == elementoNaoLinear.chaveResistiva.vLim && tempo != 0.0)
      {
        Yn[elementoNaoLinear.a][elementoNaoLinear.a]+=g_anterior[contador];
        Yn[elementoNaoLinear.b][elementoNaoLinear.b]+=g_anterior[contador];
        Yn[elementoNaoLinear.a][elementoNaoLinear.b]-=g_anterior[contador];
        Yn[elementoNaoLinear.b][elementoNaoLinear.a]-=g_anterior[contador];
      }
      else
      {
        z = 0;
        g = elementoNaoLinear.chaveResistiva.gon;
        Yn[elementoNaoLinear.a][elementoNaoLinear.a]+=g;
        Yn[elementoNaoLinear.b][elementoNaoLinear.b]+=g;
        Yn[elementoNaoLinear.a][elementoNaoLinear.b]-=g;
        Yn[elementoNaoLinear.b][elementoNaoLinear.a]-=g;
        g_anterior[contador]  = g;
        z_anterior[contador] = z;
      }
    }
    else if (elementoNaoLinear.nome[0] == 'N')
    {
      if (elementoNaoLinear.a == 0)
        newtonRaphsonAnterior[elementoNaoLinear.a] = 0;
      else if (elementoNaoLinear.b == 0)
        newtonRaphsonAnterior[elementoNaoLinear.b] = 0;

      tensaoAtual = newtonRaphsonAnterior[elementoNaoLinear.a] - newtonRaphsonAnterior[elementoNaoLinear.b];

      if (tensaoAtual < elementoNaoLinear.resistorPartes.v2)
      {
        g = (elementoNaoLinear.resistorPartes.j2 - elementoNaoLinear.resistorPartes.j1)/(elementoNaoLinear.resistorPartes.v2 - elementoNaoLinear.resistorPartes.v1);
        z = (elementoNaoLinear.resistorPartes.j2 - g*elementoNaoLinear.resistorPartes.v2);

      }
      else if (tensaoAtual == elementoNaoLinear.resistorPartes.v2 && tempo != 0.0)
      {
        g = g_anterior[contador];
        z = z_anterior[contador];
       }
      else if (tensaoAtual < elementoNaoLinear.resistorPartes.v3)
      {
        g = (elementoNaoLinear.resistorPartes.j3 - elementoNaoLinear.resistorPartes.j2)/(elementoNaoLinear.resistorPartes.v3 - elementoNaoLinear.resistorPartes.v2);
        z = (elementoNaoLinear.resistorPartes.j3 - g*elementoNaoLinear.resistorPartes.v3);
      }
      else if (tensaoAtual == elementoNaoLinear.resistorPartes.v3 && tempo != 0.0)
      {
        g = g_anterior[contador];
        z = z_anterior[contador];
      }
      else
      {
        g = (elementoNaoLinear.resistorPartes.j4 - elementoNaoLinear.resistorPartes.j3)/(elementoNaoLinear.resistorPartes.v4 - elementoNaoLinear.resistorPartes.v3);
        z = (elementoNaoLinear.resistorPartes.j4 - g*elementoNaoLinear.resistorPartes.v4);
      }
      Yn[elementoNaoLinear.a][elementoNaoLinear.a]+=g;
      Yn[elementoNaoLinear.b][elementoNaoLinear.b]+=g;
      Yn[elementoNaoLinear.a][elementoNaoLinear.b]-=g;
      Yn[elementoNaoLinear.b][elementoNaoLinear.a]-=g;
      Yn[elementoNaoLinear.a][numeroVariaveis+1]-=z;
      Yn[elementoNaoLinear.b][numeroVariaveis+1]+=z;
      g_anterior[contador]  = g;
      z_anterior[contador] = z;
    } /*resistor linear por partes*/
    #ifdef  DEBUG
      printf("Sistema remontado apos iteracao do NR\n");
      MostrarSistema();
    #endif
  }/*for*/
}/*MontarNewtonRaphson*/

unsigned TestarConvergenciaNR () /*teste de convergencia*/
{
  unsigned i;

	for (i=1; i<=numeroVariaveis;i++)
	{
		// if (fabs(Yn[i][numeroVariaveis+1]) > X_ERRO)
		// {
		// 	erros[i] = X_ERRO*fabs((Yn[i][numeroVariaveis+1]-newtonRaphsonAnterior[i])/Yn[i][numeroVariaveis+1]);
		// }
		// else
		// {
			erros[i] = fabs(Yn[i][numeroVariaveis+1]-newtonRaphsonAnterior[i]);
		//}

		if (erros[i] > MAX_ERRO_NR)
			return 0;
	}
	return 1;
}

void MontarEstampasGMin() /*acho que faz sentido, tem que testar*/
{
  //printf("Condutancia atual: %lg\n", condutanciaAtual);
  double condutancia = gminAtual;
  unsigned i, k, l, contador, no1, no2, no1Atual, no2Atual;
  double erroAtual, tensao1, tensao2;
  double fonteGMin = 0;
  unsigned vetorErros[MAX_NOS + 1];
  char tipo;
  contador = 0;
  for (i=1; i <= numeroVariaveis; i++) /*achando os nos que tem tensoes nodais com erro maior do que o permitido*/
  {
    erroAtual = erros[i];
    if (erroAtual >= MAX_ERRO_NR)
    {
      printf("Erro atual no MontarEstampas: %lg\n", erroAtual);
      vetorErros[contador] = i;
      contador += 1;
    }
  }
  for (i=0; i < contadorElementosNaoLineares; i++)
  {
    no1 = netlist[elementosNaoLineares[i]].a;
    no2 = netlist[elementosNaoLineares[i]].b;
    tipo = netlist[elementosNaoLineares[i]].nome[0];

    if (tipo == 'N') /* se resistor linear por partes coloca uma fonte de tensao em serie com a condutancia*/
    {
      tensao1 = netlist[elementosNaoLineares[i]].resistorPartes.v2;
      tensao2 = netlist[elementosNaoLineares[i]].resistorPartes.v3;
      fonteGMin = (tensao1 + tensao2)/2;
      // printf("Fonte Gmin: %lg\n", fonteGMin);
    }

    if (no1 == 0 || no2 == 0)
    {
      //printf("Entrei no if de um dos nos ser terra\n");
      for (k=0; k < contador; k++)
      {
        if (no1 == vetorErros[k])
        {
          Yn[no1][no1] += condutancia;
          Yn[no1][numeroVariaveis + 1] += fonteGMin*condutancia;
          break;
        }

        else if (no2 == vetorErros[k])
        {
          Yn[no2][no2] += condutancia;
          Yn[no2][numeroVariaveis + 1] -= fonteGMin*condutancia;
          break;
        }
      }
    }
    else
    {
      for (k=0; k < contador; k++)
      {
        no1Atual = vetorErros[k];
        for (l=k+1; l < contador; l++)
        {
          no2Atual = vetorErros[l];
          if ((no1Atual == no1 && no2Atual == no2) || (no1Atual == no2 && no2Atual == no1))
          {
            if (no1Atual == no1)
            {
              /*Monta estampa*/
              Yn[no1Atual][no1Atual] += condutancia;
              Yn[no2Atual][no2Atual] += condutancia;
              Yn[no1Atual][no2Atual] -= condutancia;
              Yn[no2Atual][no1Atual] -= condutancia;
              Yn[no1Atual][numeroVariaveis + 1] += fonteGMin*condutancia;
              Yn[no2Atual][numeroVariaveis + 1] -= fonteGMin*condutancia;
            }
            else
            {
              /*Monta estampa*/
              Yn[no1Atual][no1Atual] += condutancia;
              Yn[no2Atual][no2Atual] += condutancia;
              Yn[no1Atual][no2Atual] -= condutancia;
              Yn[no2Atual][no1Atual] -= condutancia;
              Yn[no1Atual][numeroVariaveis + 1] += fonteGMin*condutancia;
              Yn[no2Atual][numeroVariaveis + 1] -= fonteGMin*condutancia;
            }
          }
        } /*for vetorErros de dentro*/
      } /*for vetorErros de fora*/
    } /*se nenhum dos dois nos for terra*/
  } /*For elementos não lineares*/
  #ifdef  DEBUG_GMIN
    printf("Sistema remontado apos iteracao do GMin\n");
    MostrarSistema();
  #endif
} /*MontarEstampasGMin*/

unsigned ResolverNewtonRaphson (double tempo, double passo_simulacao, unsigned int pontoOperacao)
{
  unsigned convergiu;
  unsigned i,k,j, contadorFator;

  fator = 10;
  contadorFator = 1;
  ZerarSistema();
  ZerarNRAnterior();

  //Inicializacao das variaveis do sistema
  if (pontoOperacao == 1)
  {
    //Inicializacao das variaveis do sistema
  	for (i=1; i<=numeroVariaveis; i++)
  		newtonRaphsonAnterior[i] = 0.1;
  }
  else
  {
    for (i=1; i<=numeroVariaveis; i++)
      newtonRaphsonAnterior[i] = solucaoAnterior[i];
  }

  for (i=1; i <= MAX_ITERACOES; i++)
  {
    /*Monta tudo que é linear e tudo que é não-linear com base no resultado anterior*/
    MontarNewtonRaphson(tempo, passo_simulacao, pontoOperacao);
    if (ResolverSistema())
    {
      getch();
      exit(0);
    }
    convergiu = TestarConvergenciaNR();
    if (convergiu == 1)
    {
       printf("Convergiu!\n");
      //ArmazenarResultadoAnterior();
      return 1;
    }
    ArmazenarNRAnterior();
  }
    // printf("Por enquanto deu merda\n");
  //}
  //return; /*provisorio para não testar gmin*/

  ZerarSistema();
  ZerarNRAnterior();

  if (pontoOperacao == 1)
  {
    //Inicializacao das variaveis do sistema
  	for (i=1; i<=numeroVariaveis; i++)
  		newtonRaphsonAnterior[i] = 0.1;
  }
  else
  {
    for (i=1; i<=numeroVariaveis; i++)
      newtonRaphsonAnterior[i] = solucaoAnterior[i];
  }

  i = 1;
  gminAtual = GMIN_INICIAL;

  while (gminAtual > GMIN_MINIMA)
  {
    for (i=1; i <= MAX_ITERACOES; i++)
    {
      /*Monta tudo que é linear e tudo que é não-linear com base no resultado anterior*/
      MontarNewtonRaphson(tempo, passo_simulacao, pontoOperacao);
      MontarEstampasGMin();
      if (ResolverSistema())
      {
        getch();
        exit(0);
      }
      convergiu = TestarConvergenciaNR();
      if (convergiu == 1)
        break;
      ArmazenarNRAnterior();
    }
    if (convergiu)
    {
      fator = 10;
      gminAtual /= fator;
    }
    else
    {
      fator = sqrt(fator);
      gminAtual *= fator;
    }
    if (fator < GMIN_MENOR_FATOR)
    {
      printf("Nao converge nem com Gmin Stepping\n");
      return 0;
    }
  } /*while gmin*/
  return 1;

}/*ResolverNewtonRaphson*/

void CalcularMemorias (unsigned pontoOperacao, double passo_simulacao)
{
  unsigned contador;
  elemento *elementoCL;
  double resistencia, tensaoAtual, correnteAtual;
  for (contador = 0; contador < contadorElementosVariantes; contador++)
  {
    elementoCL = &(netlist[elementosVariantes[contador]]);
    if (elementoCL->nome[0] == 'C')
    {
      /*tratamento do terra*/
      if (elementoCL->a == 0)
        solucaoAnterior[elementoCL->a] = 0;
      if (elementoCL->b == 0)
        solucaoAnterior[elementoCL->b] = 0;
      /*tratamento do terra*/

      tensaoAtual = solucaoAnterior[elementoCL->a] - solucaoAnterior[elementoCL->b];

      if (pontoOperacao == 1)
      {
        correnteAtual = tensaoAtual/PO_CAPACITOR;
        //printf("Corrente capacitor: %lg\n", netlist[elementos[contador]].jt0);
      }
      else
      {
        resistencia = passo_simulacao/(2*elementoCL->valor);
        correnteAtual = tensaoAtual/(resistencia) - (1/resistencia)*(elementoCL->vt0) - elementoCL->jt0;
      }
    }/*capacitor*/

    else if (elementoCL->nome[0] == 'L') /*ainda pode dar merda*/
    {
      /*tratamento do terra*/
      if (elementoCL->a == 0)
        solucaoAnterior[elementoCL->a] = 0;
      if (elementoCL->b == 0)
        solucaoAnterior[elementoCL->b] = 0;
      /*tratamento do terra*/

      tensaoAtual = solucaoAnterior[elementoCL->a] - solucaoAnterior[elementoCL->b];

      if (pontoOperacao == 1)
      {
        correnteAtual = tensaoAtual/PO_INDUTOR;
        //printf("Corrente capacitor: %lg\n", netlist[elementos[contador]].jt0);
      }
      else
      {
        resistencia = (2*elementoCL->valor)/passo_simulacao;
        correnteAtual = tensaoAtual/resistencia + (1/resistencia)*(elementoCL->vt0) + elementoCL->jt0;
      }
    }/*indutor*/
    elementoCL->vt0 = tensaoAtual;
    elementoCL->jt0 = correnteAtual;
  }/*for*/
}/*calcular jt0 e vt0*/

/* Rotina que monta as estampas dos elementos invariantes do circuito */
void MontarEstampasInvariantes()
{
  //ZerarSistema();
	for (i=1; i<=numeroElementos; i++)
	{
		tipo=netlist[i].nome[0];
		if (tipo=='R')
		{
			g = 1/netlist[i].valor;
			YnInvariantes[netlist[i].a][netlist[i].a] += g;
			YnInvariantes[netlist[i].b][netlist[i].b] += g;
			YnInvariantes[netlist[i].a][netlist[i].b] -= g;
			YnInvariantes[netlist[i].b][netlist[i].a] -= g;
		}
		else if (tipo=='G')
		{
			g = netlist[i].valor;
			YnInvariantes[netlist[i].a][netlist[i].c] += g;
			YnInvariantes[netlist[i].b][netlist[i].d] += g;
			YnInvariantes[netlist[i].a][netlist[i].d] -= g;
			YnInvariantes[netlist[i].b][netlist[i].c] -= g;
		}

    /*So preenche se for DC*/
    else if (strcmp(netlist[i].tipo_fonte, "DC") == 0)
    {
      netlist[i].valor = netlist[i].fonte_dc.valor;
      g = netlist[i].valor;
      if (tipo =='I')
      {
        YnInvariantes[netlist[i].a][numeroVariaveis+1]-=g;
        YnInvariantes[netlist[i].b][numeroVariaveis+1]+=g;
      }
      else if (tipo == 'V')
      {
        YnInvariantes[netlist[i].a][netlist[i].x]+=1;
        YnInvariantes[netlist[i].b][netlist[i].x]-=1;
        YnInvariantes[netlist[i].x][netlist[i].a]-=1;
        YnInvariantes[netlist[i].x][netlist[i].b]+=1;
        YnInvariantes[netlist[i].x][numeroVariaveis+1]-=g;
      }
    }
		else if (tipo=='E')
		{
			g = netlist[i].valor;
			YnInvariantes[netlist[i].a][netlist[i].x] += 1;
			YnInvariantes[netlist[i].b][netlist[i].x] -= 1;
			YnInvariantes[netlist[i].x][netlist[i].a] -= 1;
			YnInvariantes[netlist[i].x][netlist[i].b] += 1;
			YnInvariantes[netlist[i].x][netlist[i].c] += g;
			YnInvariantes[netlist[i].x][netlist[i].d] -= g;
		}
		else if (tipo=='F')
		{
			g = netlist[i].valor;
			YnInvariantes[netlist[i].a][netlist[i].x] += g;
			YnInvariantes[netlist[i].b][netlist[i].x] -= g;
			YnInvariantes[netlist[i].c][netlist[i].x] += 1;
			YnInvariantes[netlist[i].d][netlist[i].x] -= 1;
			YnInvariantes[netlist[i].x][netlist[i].c] -= 1;
			YnInvariantes[netlist[i].x][netlist[i].d] += 1;
		}
		else if (tipo=='H')
		{
			g = netlist[i].valor;
			YnInvariantes[netlist[i].a][netlist[i].y] += 1;
			YnInvariantes[netlist[i].b][netlist[i].y] -= 1;
			YnInvariantes[netlist[i].c][netlist[i].x] += 1;
			YnInvariantes[netlist[i].d][netlist[i].x] -= 1;
			YnInvariantes[netlist[i].y][netlist[i].a] -= 1;
			YnInvariantes[netlist[i].y][netlist[i].b] += 1;
			YnInvariantes[netlist[i].x][netlist[i].c] -= 1;
			YnInvariantes[netlist[i].x][netlist[i].d] += 1;
			YnInvariantes[netlist[i].y][netlist[i].x] += g;
		}
		else if (tipo=='O')
		{
			YnInvariantes[netlist[i].a][netlist[i].x] += 1;
			YnInvariantes[netlist[i].b][netlist[i].x] -= 1;
			YnInvariantes[netlist[i].x][netlist[i].c] += 1;
			YnInvariantes[netlist[i].x][netlist[i].d] -= 1;
		}
    else if (tipo == 'K')
    {
      printf("Numero de espiras: %lg\n", netlist[i].numeroEspiras);
      YnInvariantes[netlist[i].a][netlist[i].x] -= netlist[i].numeroEspiras;
			YnInvariantes[netlist[i].b][netlist[i].x] += netlist[i].numeroEspiras;
      YnInvariantes[netlist[i].c][netlist[i].x] += 1;
      YnInvariantes[netlist[i].d][netlist[i].x] -= 1;
			YnInvariantes[netlist[i].x][netlist[i].a] += netlist[i].numeroEspiras;
			YnInvariantes[netlist[i].x][netlist[i].b] -= netlist[i].numeroEspiras;
      YnInvariantes[netlist[i].x][netlist[i].c] -= 1;
			YnInvariantes[netlist[i].x][netlist[i].d] += 1;
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
	}/*for*/
}/*MontarEstampasInvariantes*/

void ListarTudo ()
{
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
}/*ListarTudo*/

void AcrescentarVariaveis()
{
  /* Acrescenta variaveis de corrente acima dos nos, anotando no netlist */
  numeroNos=numeroVariaveis;
  for (i=1; i<=numeroElementos; i++)
  {
    tipo=netlist[i].nome[0];
    if (tipo=='V' || tipo=='E' || tipo=='F' || tipo=='O' || tipo=='K')
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
}/*AcrescentarVariaveis*/

void ResolverPontoOperacao (double passo_simulacao)
{
  MontarEstampasVariantes(0, passo_simulacao, 1); /*ponto de operacao*/
  #ifdef DEBUG
    MostrarSistema();
  #endif
  if (contadorElementosNaoLineares == 0) /*se nao tem elementos nao lineares, resolve normalmente*/
  {
    printf("Entra no if de nao tem elementos nao lineares\n");
    if (ResolverSistema()) /*calculo do ponto de operacao*/
    {
      getch();
      exit(0);
    }
    ArmazenarResultadoAnterior(); /*armazeno as solucoes anteriores num vetor solucaoAnterior*/
    CalcularMemorias(1, passo_simulacao); /*armazeno as correntes nos capacitores e as tensoes nos indutores do resultado anterior*/
  }
  else /*se tiver elementos nao lineares, resolve com NR*/
  {
    if (ResolverNewtonRaphson(0, passo_simulacao, 1))
      ArmazenarResultadoAnterior();
    /*Newton-Raphson que depende da solucao anterior*/
  }
} /*ResolverPontoOperacao*/

int main(void)
{
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
  numeroElementos=0;
  numeroVariaveis=0;
  double tempo_atual = 0.0;
  unsigned convergencia = 0;
  unsigned contadorPasso = 1;
  strcpy(lista[0],"0");
  printf("Nome do arquivo com o netlist (ex: mna.net): ");
  scanf("%50s",nomeArquivo);
  arquivo=fopen(nomeArquivo,"r");
  if (arquivo==0)
  {
    printf("Arquivo %s inexistente\n",nomeArquivo);
    goto denovo;
  }

  LerNetlist(arquivo);
  AcrescentarVariaveis();
  ListarTudo();

  FILE *arquivoSaida = fopen(NOME_ARQUIVO_SAIDA, "w");
  /*Escreve o header do arquivo de saida*/
  fprintf(arquivoSaida, "%s", "t");
  for (i=1; i<=numeroVariaveis; i++)
    fprintf(arquivoSaida, " %s", lista[i]);
  fprintf(arquivoSaida, "\n");

  if (passosPorPonto != 1)
  {
    // printf("Passos por ponto %u\n", passosPorPonto);
    passo_simulacao = passo_simulacao/passosPorPonto;
    // printf("Passo atual: %lg\n", passoTeste);

  }
  ZerarSistema();
  MontarEstampasInvariantes();
  CopiarEstampaInvariante();
  printf("Sistema apos estampas invariantes:\n");
  MostrarSistema();
  if (temCapacitorOuIndutor == 1 || contadorElementosNaoLineares != 0)
  {
    ResolverPontoOperacao(passo_simulacao);
    MostrarSolucaoAtual();
  }

  fprintf(arquivoSaida,"%lg", tempo_atual);
  for (i=1; i<=numeroVariaveis; i++)
  {
    if (contadorElementosNaoLineares != 0)
      fprintf(arquivoSaida," %lg", solucaoAnterior[i]);
    else
      fprintf(arquivoSaida," %lg", Yn[i][numeroVariaveis+1]);
  }
  fprintf(arquivoSaida,"\n");


  /*Analise no tempo*/
  for (tempo_atual = passo_simulacao; tempo_atual <= tempo_simulacao; tempo_atual += passo_simulacao) /*começa em 0 ou em 0 + passo?*/
  {

    if (contadorElementosNaoLineares != 0)
    {
      /*Newton-Raphson com parametros do sistema inicial (ponto de operacao)*/
      convergencia = ResolverNewtonRaphson(tempo_atual, passo_simulacao, 0);
      if (convergencia)
        ArmazenarResultadoAnterior();

      if (temCapacitorOuIndutor == 1) /*COMO FAZER ISSO DE UMA FORMA MAIS ESPERTA?*/
        CalcularMemorias(0, passo_simulacao);  /*armazeno as correntes nos capacitores e as tensoes nos indutores do resultado anterior*/
        //ZerarResultadoAnterior();
    }
    else
    {
      MontarEstampasVariantes(tempo_atual, passo_simulacao, 0);
      /* Resolve o sistema */
      if (ResolverSistema())
      {
        getch();
        exit(0);
      }
      if (temCapacitorOuIndutor == 1) /*COMO FAZER ISSO DE UMA FORMA MAIS ESPERTA?*/
      {
        ArmazenarResultadoAnterior();
        CalcularMemorias(0, passo_simulacao);  /*armazeno as correntes nos capacitores e as tensoes nos indutores do resultado anterior*/
      }
    }

    if (contadorPasso == passosPorPonto || passosPorPonto == 1)
    {
      /*Escreve no arquivo de saida*/
      fprintf(arquivoSaida,"%lg", tempo_atual);
      for (i=1; i<=numeroVariaveis; i++)
      {
        if (contadorElementosNaoLineares != 0)
          fprintf(arquivoSaida," %lg", solucaoAnterior[i]);
        else
          fprintf(arquivoSaida," %lg", Yn[i][numeroVariaveis+1]);
      }
      fprintf(arquivoSaida,"\n");
      contadorPasso = 1;
    }

    contadorPasso++;

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
