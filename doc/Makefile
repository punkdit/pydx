
# This uses docutils (rst2html) and epydoc version 3.0
# http://docutils.sourceforge.net/
# http://epydoc.sourceforge.net/

all: api-doc # guide 

guide:
	rst2html.py --stylesheet-path=./style.css ./getting-started.txt ./getting-started.html 
api-doc:
	epydoc -v -o api --name pydx --docformat=epytext \
		pydx.scalar pydx.scalar.symbolic   pydx.scalar.mpfi pydx.scalar.fmath\
		pydx.mjet pydx.tensor pydx.field pydx.transform pydx.ode pydx.geodesic \
		pydx.metric 

