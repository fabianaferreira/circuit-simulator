package system

import (
	"fmt"
	"math"
)

const tolg = 1e-9

type system struct {
	elements, variables, nodes int
	errors                     []float64
	matrix                     [][]float64
	previousResult             []float64
	previousNewtonRaphson      []float64
	newtonRaphsonMatrix        [][]float64
	invariantPrints            [][]float64
}

func NewSystemSolver(elems, vars, nodes int) *system {
	return &system{
		elements:              elems,
		variables:             vars,
		nodes:                 nodes,
		errors:                make([]float64, vars),
		matrix:                make([][]float64, vars+1),
		newtonRaphsonMatrix:   make([][]float64, vars+1),
		previousResult:        make([]float64, vars+2),
		previousNewtonRaphson: make([]float64, vars+2),
		invariantPrints:       make([][]float64, vars+1),
	}
}

func (s *system) SolveSystem() error {
	var i, j, l, a int
	var t, p float64

	for i = 1; i <= s.variables; i++ {
		t = 0.0
		a = i

		for l = i; l <= s.variables; l++ {

			if math.Abs(s.matrix[l][i]) > math.Abs(t) {
				a = l
				t = s.matrix[l][i]
			}
		}

		if i != a {
			for l := 1; l <= s.variables+1; l++ {
				p = s.matrix[i][l]
				s.matrix[i][l] = s.matrix[a][l]
				s.matrix[a][l] = p
			}
		}
		if math.Abs(t) < tolg {
			fmt.Println("Irregular system")

			panic("Irregular system")
		}

		for j = s.variables + 1; j > 0; j-- {
			s.matrix[i][j] /= t
			p = s.matrix[i][j]
			if p != 0 {
				for l := 1; l <= s.variables; l++ {
					if l != i {
						s.matrix[l][j] -= s.matrix[l][i] * p
					}
				}
			}
		}
	}
	return nil
}

func (s *system) SavePreviousResult() {
	for i := 1; i <= s.variables; i++ {
		s.previousResult[i] = s.matrix[i][s.variables+1] // last column from matrix: system's solution
	}
	s.previousResult[0] = 0
}

func (s *system) InitializeErrorVector() {
	for index := range s.errors {
		s.errors[index] = 1
	}
}

func (s *system) ResetErrorVector() {
	for index := range s.errors {
		s.errors[index] = 0
	}
}
