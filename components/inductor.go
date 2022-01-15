package components

const inductorP0 = 1e-9

func (e *Element) CalculateInductorMemory(operationPoint uint, simulationStep, previousA, previousB float64) {
	var tension, current, resistance float64
	tension = previousA - previousB
	if operationPoint == 1 {
		current = tension / inductorP0
	} else {
		resistance = 2 * e.Value / simulationStep
		current = tension/resistance + (1/resistance)*e.Vt0 + e.Jt0
	}

	e.Vt0 = tension
	e.Jt0 = current
}
