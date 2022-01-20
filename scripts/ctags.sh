#!/bin/sh
#
# Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>
# Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg
# modified on 2022-01-20 by stefan.kurtz@uni-hamburg.de
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

# try to differntiate between Exuberant Ctags and BSD Ctags
ctags --version > /dev/null

if [ $? -eq 0 ]
then
  exuberantopts="--c++-kinds=+p --fields=+iaSKlm --extra=+q"
fi

# make a new tags file
ctags -w ${exuberantopts} src/*/*.[ch]pp \
                          testsuite/*.cpp
