package main

// #define versao              "1.0"
// #define MAX_LINHA           80
// #define MAX_TIPO_FONTE      5
// #define MAX_NOME            11
// #define MAX_ELEM            50
// #define MAX_NOS             50
// #define TOLG                1e-9
// #define PI                  acos(-1.0)
// #define PO_CAPACITOR        1e9
// #define PO_INDUTOR          1e-9
// #define MAX_ERRO_NR         1e-6
// #define X_ERRO              1
// #define MAX_ITERACOES       100
// #define MAX_INICIALIZACOES  10
// #define GMIN_INICIAL        1.1
// #define GMIN_MINIMA         1e-12
// #define GMIN_MENOR_FATOR    1.1
// #define INFINITO            1e12
// #define NOME_ARQUIVO_SAIDA "saida_simulacao.tab"

const (
	tolg               = 1e-9
	capacitorOP        = 1e9
	inductorOP         = 1e-9
	newRaphsonMaxError = 1e-6
	xError             = 1
	maxIterations      = 100
	initialGMin        = 1e-12
	minorFactorGmin    = 1.1
	outputFilename     = "simulation_output.tab"
)
