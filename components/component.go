package components

/*Elemento possui atributos de fontes, pois há tipos diferentes, com parametros diferentes*/
type Element struct /* Elemento do netlist */
{
	Name             string
	Value            float64
	FontType         string
	A, B, C, D, X, Y int
	Sine             Sine
	DcSource         Dc
	PulseSource      Pulse
	ResistiveKey     Key
	Resistor         Resistor
	SpiresQnt        uint64
	Jt0, Vt0         float64
}
