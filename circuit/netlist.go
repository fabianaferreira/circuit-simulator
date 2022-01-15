package circuit

import (
	"bufio"
	"fmt"
	"os"
	"strings"

	"github.com/FabianaFerreira/circuit-simulator/components"
	"github.com/FabianaFerreira/circuit-simulator/utils"
)

const (
	maxElements   = 50
	maxNodes      = 50
	maxLines      = 80
	maxFonts      = 5
	maxNameLength = 11
)

func getNetlistLines(filename string) ([]string, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("failed to open file: %s", err)
	}

	defer f.Close()

	fileScanner := bufio.NewScanner(f)
	fileScanner.Split(bufio.ScanLines)

	var fileTextLines []string
	for fileScanner.Scan() {
		fileTextLines = append(fileTextLines, fileScanner.Text())
	}

	return fileTextLines, nil
}

func ReadNetlist(filename string) error {
	lines, err := getNetlistLines(filename)
	if err != nil {
		return fmt.Errorf("error while reading netlist: %s", err)
	}

	c := NewCircuit()
	fmt.Printf("Netlist title: %s", lines[0])
	// First line is the title
	lines = lines[1:]

	for _, line := range lines {
		c.elementsQuantity++
		if c.elementsQuantity > maxElements {
			return fmt.Errorf("maximum number of elements is: %v", maxElements)
		}

		elemName, elemType, err := getElementNameAndType(line)
		if err != nil {
			return err
		}

		// Removes element name
		// line =

		c.netlist = append(c.netlist, components.Element{Name: elemName})

		elem, err := processLine(elemType, line[1:], c)
		if err != nil {
			return err
		}

		c.netlist[c.elementsQuantity] = *elem
	}

	return nil
}

func getElementNameAndType(str string) (string, string, error) {
	substrings := strings.Fields(str)
	elemName := substrings[0]
	if len(elemName) > maxNameLength {
		return "", "", fmt.Errorf("error while reading netlist: element name size (%v) is greater than allowed (%v)", len(elemName), maxNameLength)
	}

	return substrings[0], string(substrings[0][0]), nil
}

func processLine(elemType, line string, c *Circuit) (*components.Element, error) {
	var na, nb, nc, nd, name string
	currentElem := &components.Element{}

	switch elemType {
	case "R":
		if _, err := fmt.Sscanf(line, "%s%10s%10s%f", &name, &na, &nb, &currentElem.Value); err != nil {
			return nil, err
		}
		nodes, err := c.identifyNodes(na, nb)
		if err != nil {
			return nil, err
		}

		utils.Unpack(nodes, &currentElem.A, &currentElem.B)

	case "C":
		fmt.Sscanf(line, "%s%10s%10s%f", &name, &na, &nb, &currentElem.Value)
		nodes, err := c.identifyNodes(na, nb)
		if err != nil {
			return nil, err
		}

		utils.Unpack(nodes, &currentElem.A, &currentElem.B)

		c.hasComponentWithMemory = true
		c.variantElements = append(c.variantElements, c.elementsQuantity)

	case "L":
		fmt.Sscanf(line, "%s%10s%10s%f", &name, &na, &nb, &currentElem.Value)
		nodes, err := c.identifyNodes(na, nb)
		if err != nil {
			return nil, err
		}

		utils.Unpack(nodes, &currentElem.A, &currentElem.B)

		c.hasComponentWithMemory = true
		c.variantElements = append(c.variantElements, c.elementsQuantity)

	case "I", "V":
		fmt.Sscanf(line, "%10s%10s%5s", &na, &nb, &currentElem.FontType)
		nodes, err := c.identifyNodes(na, nb)
		if err != nil {
			return nil, err
		}

		utils.Unpack(nodes, &currentElem.A, &currentElem.B)

		switch fontType := currentElem.FontType; fontType {
		case "DC":
			fmt.Sscanf(line, "%*10s%*10s%*5s%lg", currentElem.DcSource.Value)

		case "SIN":
			fmt.Sscanf(line, "%*10s%*10s%*5s%lg%lg%lg%lg%lg%lg%i",
				currentElem.Sine.Amplitude,
				currentElem.Sine.Freq,
				currentElem.Sine.Delay,
				currentElem.Sine.Damping,
				currentElem.Sine.Phase,
				currentElem.Sine.Cycles,
			)
			c.variantElements = append(c.variantElements, c.elementsQuantity)

		case "PULSE":
			fmt.Sscanf(line, "%*10s%*10s%*5s%lg%lg%lg%lg%lg%lg%lg%i",
				currentElem.PulseSource.Amplitude1,
				currentElem.PulseSource.Amplitude2,
				currentElem.PulseSource.Delay,
				currentElem.PulseSource.RisingTime,
				currentElem.PulseSource.FallingTime,
				currentElem.PulseSource.Duration,
				currentElem.PulseSource.Period,
				currentElem.PulseSource.Cycles,
			)
			c.variantElements = append(c.variantElements, c.elementsQuantity)
		}

	case "G", "E", "F", "H":
		fmt.Sscanf(line, "%10s%10s%10s%10s%lg", na, nb, nc, nd, currentElem.Value)
		nodes, err := c.identifyNodes(na, nb, nc, nd)
		if err != nil {
			return nil, err
		}

		utils.Unpack(nodes,
			&currentElem.A,
			&currentElem.B,
			&currentElem.C,
			&currentElem.D,
		)

	case "O":
		fmt.Sscanf(line, "%10s%10s%10s%10s", na, nb, nc, nd)
		nodes, err := c.identifyNodes(na, nb, nc, nd)
		if err != nil {
			return nil, err
		}

		utils.Unpack(nodes,
			&currentElem.A,
			&currentElem.B,
			&currentElem.C,
			&currentElem.D,
		)

	case "K":
		fmt.Sscanf(line, "%10s%10s%10s%10s%lg", na, nb, nc, nd, currentElem.SpiresQnt)
		nodes, err := c.identifyNodes(na, nb, nc, nd)
		if err != nil {
			return nil, err
		}

		utils.Unpack(nodes,
			&currentElem.A,
			&currentElem.B,
			&currentElem.C,
			&currentElem.D,
		)

	case "N":
		fmt.Sscanf(line, "%10s%10s%lg%lg%lg%lg%lg%lg%lg%lg", na, nb,
			&currentElem.Resistor.V1,
			&currentElem.Resistor.J1,
			&currentElem.Resistor.V2,
			&currentElem.Resistor.J2,
			&currentElem.Resistor.V3,
			&currentElem.Resistor.J3,
			&currentElem.Resistor.V4,
			&currentElem.Resistor.J4,
		)

		nodes, err := c.identifyNodes(na, nb, nc, nd)
		if err != nil {
			return nil, err
		}

		utils.Unpack(nodes, &currentElem.A, &currentElem.B)
		c.nonLinearElements = append(c.nonLinearElements, c.elementsQuantity)

	case "$":
		fmt.Sscanf(line, "%10s%10s%10s%10s%lg%lg%lg", na, nb, nc, nd,
			currentElem.ResistiveKey.Gon,
			currentElem.ResistiveKey.Goff,
			currentElem.ResistiveKey.VLim)

		nodes, err := c.identifyNodes(na, nb, nc, nd)
		if err != nil {
			return nil, err
		}

		utils.Unpack(nodes, &currentElem.A, &currentElem.B)
		c.nonLinearElements = append(c.nonLinearElements, c.elementsQuantity)

	case "*":
		fmt.Printf("Comentario: %s \n", line)
		c.elementsQuantity--
	case ".":
		if strings.Compare(currentElem.Name, ".TRAN") == 0 {
			fmt.Sscanf(line, "%lg%lg%*10s%u", c.simulationTime, c.simulationSteps, c.stepsPerPoint)
			c.elementsQuantity--
		}
	default:
		fmt.Printf("Elemento desconhecido: %s\n", line)
		os.Exit(0)
	}

	return currentElem, nil
}
