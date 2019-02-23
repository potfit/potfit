CXX = g++
CXXFLAGS = $(shell pkg-config libkim-api-v2 --cflags)
LDFLAGS = $(shell pkg-config libkim-api-v2 --libs)

default: kim_read_model_props

clean:
	rm -rf kim_read_model_props

kim_read_model_props: kim_read_model_props.cpp
	@echo "Building kim_read_model_props tool"
	@$(CXX) $(CXXFLAGS) -fsanitize=address -std=c++17 kim_read_model_props.cpp $(LDFLAGS) -o kim_read_model_props
