# run project locally
run:
	@go run main.go
	
# run tests for the project
test:
	@go test ./... -coverprofile cover.out -v -race

coverage: test
	@go tool cover -html=cover.out