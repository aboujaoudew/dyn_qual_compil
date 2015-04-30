MAINREP=tool/
MODELREP=$(MAINREP)input/
EXAMPLES=$(notdir $(wildcard $(MODELREP)*))
OUTPUTREP=$(MAINREP)result/
DOT=$(wildcard $(OUTPUTREP)*/*.dot)
PDF=$(DOT:%.dot=%.pdf)

DOTENGINE=dot
MATLAB=matlab

all: $(EXAMPLES)
	make -C . pdf 


pdf: $(PDF)

%: $(MODELREP)%
	cd $(MAINREP) ; rm -rf tmp.m ; echo "frontend('$@')" > tmp.m ; $(MATLAB) tmp.m 

%.pdf: %.dot
	$(DOTENGINE) -Tpdf $< -o $@
