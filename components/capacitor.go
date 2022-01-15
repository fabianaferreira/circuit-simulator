package components

const capacitorP0 = 1e9

func (e *Element) CalculateCapacitorMemory(operationPoint uint, simulationStep, previousA, previousB float64) {
	var tension, current, resistance float64
	tension = previousA - previousB
	if operationPoint == 1 {
		current = tension / capacitorP0
	} else {
		resistance = simulationStep / (2 * e.Value)
		current = tension/resistance - (1/resistance)*e.Vt0 - e.Jt0
	}

	e.Vt0 = tension
	e.Jt0 = current
}
