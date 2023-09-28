################################################################################
#                                   LICENSE                                    #
################################################################################
#   This file is part of rss_ringoccs.                                         #
#                                                                              #
#   rss_ringoccs is free software: you can redistribute it and/or modify       #
#   it under the terms of the GNU General Public License as published by       #
#   the Free Software Foundation, either version 3 of the License, or          #
#   (at your option) any later version.                                        #
#                                                                              #
#   rss_ringoccs is distributed in the hope that it will be useful,            #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#   GNU General Public License for more details.                               #
#                                                                              #
#   You should have received a copy of the GNU General Public License          #
#   along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.     #
################################################################################
#   Author:     Ryan Maguire                                                   #
#   Date:       Sepemtber 15, 2022                                             #
################################################################################

TARGET_LIB := librssringoccs.so
BUILD_DIR := ./build
SRC_DIRS := ./src

ifdef OMP
CFLAGS := $(EXTRA_FLAGS) -fopenmp -I../ -O3 -fPIC -flto -DNDEBUG -c
LFLAGS := $(EXTRA_LFLAGS) -fopenmp -O3 -flto -shared -lm -ltmpl
else
CFLAGS := $(EXTRA_FLAGS) -I../ -O3 -fPIC -flto -DNDEBUG -c
LFLAGS := $(EXTRA_LFLAGS) -O3 -flto -shared -lm -ltmpl
endif

CWARN := -Wall -Wextra -Wpedantic
SRCS := $(shell find $(SRC_DIRS) -name "*.c")
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

.PHONY: clean install uninstall all

all: $(BUILD_DIR) $(TARGET_LIB)

$(TARGET_LIB): $(OBJS)
	$(CC) $(OBJS) $(LFLAGS) -o $@

$(BUILD_DIR)/%.c.o: %.c
	$(CC) $(CWARN) $(CFLAGS) $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)/src/csv_tools/
	mkdir -p $(BUILD_DIR)/src/fresnel_kernel/
	mkdir -p $(BUILD_DIR)/src/fresnel_transform/
	mkdir -p $(BUILD_DIR)/src/history/
	mkdir -p $(BUILD_DIR)/src/occultation_geometry/
	mkdir -p $(BUILD_DIR)/src/reconstruction/
	mkdir -p $(BUILD_DIR)/src/tau/

clean:
	rm -rf $(BUILD_DIR)
	rm -f $(TARGET_LIB)

install:
	mkdir -p /usr/local/lib/
	mkdir -p /usr/local/include/rss_ringoccs/
	cp -r ./include /usr/local/include/rss_ringoccs/
	cp $(TARGET_LIB) /usr/local/lib/$(TARGET_LIB)

uninstall:
	rm -rf $(BUILD_DIR)
	rm -f $(TARGET_LIB)
	rm -rf /usr/local/include/rss_ringoccs/
	rm -f /usr/local/lib/$(TARGET_LIB)

-include $(DEPS)

