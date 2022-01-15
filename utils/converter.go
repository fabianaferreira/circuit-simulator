package utils

import "strconv"

func StringToFloat64(str string) (float64, error) {
	return strconv.ParseFloat(str, 64)
}
