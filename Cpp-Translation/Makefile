CXX      := clang++
CXXFLAGS := -std=c++23 -g -Wall -Wextra -O3 -march=native -MMD -MP -Iinclude
LDFLAGS  :=

BUILD_DIR := build
BIN       := fmm-run

SRC := $(shell find src -name '*.cpp')
OBJ := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(SRC))
DEP := $(OBJ:.o=.d)

.PHONY: all debug run clean

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@ $(LDFLAGS)

$(BUILD_DIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rebuilds everything with sanitizers + debug symbols, no optimization.
debug: CXXFLAGS := -std=c++23 -Wall -Wextra -g -O0 -fsanitize=address,undefined -Iinclude
debug: clean $(BIN)

run: $(BIN)
	./$(BIN)

clean:
	rm -rf $(BUILD_DIR) $(BIN)

-include $(DEP)
