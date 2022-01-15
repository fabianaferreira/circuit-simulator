package utils

func Unpack(s []int, vars ...*int) {
	for i, value := range s {
		*vars[i] = value
	}
}
