package system

func (s *system) CopyInvariantPrints() {
	for i := 0; i <= s.variables; i++ {
		for j := 0; j <= s.variables+1; j++ {
			s.matrix[i][j] = s.invariantPrints[i][j]
		}
	}
}

func (s *system) SaveNewtonRaphsonPrints() /*vai armazenar a estampa do NR para adicionar o Gmin*/ {
	for i := 0; i <= s.variables; i++ {
		for j := 0; j < s.variables+1; j++ {
			s.newtonRaphsonMatrix[i][j] = s.matrix[i][j]
		}
	}
}

func (s *system) CopyNewtonRaphsonPrints() {
	for i := 0; i <= s.variables; i++ {
		for j := 0; j < s.variables+1; j++ {
			s.matrix[i][j] = s.newtonRaphsonMatrix[i][j]
		}
	}
}
