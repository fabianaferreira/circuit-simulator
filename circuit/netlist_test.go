package circuit

import (
	"testing"

	"github.com/FabianaFerreira/circuit-simulator/components"
	. "github.com/smartystreets/goconvey/convey"
)

func TestParsingLine(t *testing.T) {
	Convey("Given a line", t, func() {
		line := "R0100 1 0 1"
		elemType := "R"
		circuit := NewCircuit()

		Convey("When I try to process the line", func() {
			elem, err := processLine(elemType, line, circuit)
			So(err, ShouldBeNil)

			Convey("And the elem read from the line should match", func() {
				expectedElem := components.Element{
					Value: 1,
					A:     1,
					B:     0,
				}
				So(*elem, ShouldResemble, expectedElem)
			})
		})
	})
}
