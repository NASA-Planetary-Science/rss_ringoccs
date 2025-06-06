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
TARGET_LIB_SHARED := librssringoccs.so
TARGET_LIB_STATIC := librssringoccs.a
ifdef BUILD_STATIC
TARGET_LIB := $(TARGET_LIB_STATIC)
AR ?= ar
else
TARGET_LIB := $(TARGET_LIB_SHARED)
endif

BUILD_DIR := ./build
SRC_DIRS := ./src
LIBTMPL := ./libtmpl/libtmpl.a

CWARN := -Wall -Wextra -Wpedantic
CFLAGS := $(EXTRA_FLAGS) -I../ -I./ -O3 -fPIC -flto -DNDEBUG -c
LFLAGS := $(EXTRA_LFLAGS) -L./libtmpl/ -O3 -flto -shared
LFLAGS += -lm -l:libtmpl.a -l:libcspice.a -l:libcsupport.a

ifdef OMP
CFLAGS += -fopenmp
LFLAGS += -fopenmp
endif

SRCS := $(shell find $(SRC_DIRS) -name "*.c")
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

.PHONY: clean install uninstall all libtmpl

all: $(TARGET_LIB) $(LIBTMPL)

$(TARGET_LIB_SHARED): $(OBJS) $(LIBTMPL)
	@echo "Building librssringoccs.so ..."
	@-$(CC) $(OBJS) $(LFLAGS) -o $@

$(TARGET_LIB_STATIC): $(OBJS) $(LIBTMPL)
	@echo "Building librssringoccs.a ..."
	@-$(AR) rcs $@ $(OBJS)

$(BUILD_DIR)/%.c.o: %.c $(LIBTMPL)
	@mkdir -p $(@D)
	$(CC) $(CWARN) $(CFLAGS) $< -o $@

$(LIBTMPL):
	$(MAKE) -C libtmpl/ BUILD_STATIC=1

clean:
	rm -rf $(BUILD_DIR)
	rm -f $(TARGET_LIB)
	$(MAKE) -C libtmpl/ clean BUILD_STATIC=1

install:
	mkdir -p /usr/local/lib/
	mkdir -p /usr/local/include/rss_ringoccs/
	cp -r ./include /usr/local/include/rss_ringoccs/
	cp $(TARGET_LIB) /usr/local/lib/$(TARGET_LIB)
	mkdir -p /usr/local/lib/
	mkdir -p /usr/local/include/libtmpl/
	cp -r ./libtmpl/include /usr/local/include/libtmpl/

uninstall:
	rm -rf $(BUILD_DIR)
	rm -f $(TARGET_LIB)
	rm -rf /usr/local/include/rss_ringoccs/
	rm -rf /usr/local/include/libtmpl/
	rm -f /usr/local/lib/$(TARGET_LIB)
	$(MAKE) -C libtmpl/ uninstall
