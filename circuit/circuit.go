package circuit

import (
	"fmt"
	"strings"

	"github.com/FabianaFerreira/circuit-simulator/components"
	"github.com/FabianaFerreira/circuit-simulator/utils"
)

var (
	list = make([]string, 0, maxNodes+1)
)

type Circuit struct {
	elementsQuantity, variables, nodes int
	netlist                            []components.Element
	hasComponentWithMemory             bool
	variantElements                    []int
	nonLinearElements                  []int
	simulationTime                     int64
	simulationSteps                    int64
	stepsPerPoint                      int
}

func NewCircuit() *Circuit {
	list = append(list, "0")

	return &Circuit{
		elementsQuantity:       0,
		variables:              0,
		nodes:                  0,
		netlist:                make([]components.Element, 0, maxElements),
		hasComponentWithMemory: false,
		variantElements:        make([]int, 0),
		nonLinearElements:      make([]int, 0),
	}
}

func (c *Circuit) AddVariables() error {
	c.nodes = c.variables
	var elemType string

	for i := 1; i <= c.elementsQuantity; i++ {
		elemType = string(c.netlist[i].Name[0])

		switch t := elemType; t {
		case "H":
			c.variables += 2
			if c.variables > maxNodes {
				return fmt.Errorf("currents exceeded the number of variables allowed (%d)", maxNodes)
			}
			sx := []string{"jx", c.netlist[i].Name}
			list = append(list, strings.Join(sx, ""))
			c.netlist[i].X = c.variables - 1

			sy := []string{"jy", c.netlist[i].Name}
			list = append(list, strings.Join(sy, ""))
			c.netlist[i].Y = c.variables

		default:
			c.variables++
			if c.variables > maxNodes {
				return fmt.Errorf("currents exceeded the number of variables allowed (%d)", maxNodes)
			}
			s := []string{"j", c.netlist[i].Name}
			list = append(list, strings.Join(s, ""))
			c.netlist[i].X = c.variables
		}
	}

	return nil
}

func (c *Circuit) identifyNodes(names ...string) ([]int, error) {
	nodes := make([]int, 0)

	for _, name := range names {
		found, index := utils.Find(list, name)

		if !found {
			if c.variables == maxNodes {
				return nil, fmt.Errorf("maximum number of nodes is %v", maxNodes)
			}
			c.variables++
			list = append(list, name)
			nodes = append(nodes, c.variables)
		} else {
			nodes = append(nodes, index)
		}
	}

	return nodes, nil
}
