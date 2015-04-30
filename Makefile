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


models: $(EXAMPLES)
pdf: $(PDF)

%: $(MODELREP)%
	cd $(MAINREP) ; rm -rf tmp.m ; echo "frontend('$@')" > tmp.m ; $(MATLAB) tmp.m 

%.pdf: %.dot
	$(DOTENGINE) -Tpdf $< -o $@

help:
	@echo make: compile all the models and generate the corresponfing pdf
	@echo make models: compile all models
	@echo make pdf: convert each dot file into a pdf file 
	@echo make prozone: compile the prozone model 
	@echo make rescaling: compile the model with several time scales
