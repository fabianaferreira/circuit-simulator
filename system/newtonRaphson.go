package system

func (s *system) SavePreviousNewtonRaphson() {
	for i := 1; i <= s.variables; i++ {
		s.previousNewtonRaphson[i] = s.matrix[i][s.variables+1] // last column from matrix: system's solution
	}
}

func (s *system) ResetPreviousNewtonRaphson() {
	for i := 1; i <= s.variables; i++ {
		s.previousNewtonRaphson[i] = 0
	}
}
