package prints

import "github.com/FabianaFerreira/circuit-simulator/components"

// TODO: evaluate the best option for not copying the matrix for every component method
func AddResistor(matrix [][]float64, element *components.Element) *[][]float64 {
	g := 1 / element.Value
	matrix[element.A][element.A] += g
	matrix[element.B][element.B] += g
	matrix[element.A][element.B] -= g
	matrix[element.B][element.A] -= g

	return &matrix
}
