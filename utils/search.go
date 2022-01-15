package utils

func Find(list []string, desired string) (bool, int) {
	for index, elem := range list {
		if desired == elem {
			return true, index
		}
	}
	return false, len(list)
}
