MAINREP=main/
MODELREP=case_studies/
EXAMPLES=$(notdir $(wildcard $(MODELREP)*))
OUTPUTREP=$(wildcard $(MODELREP)*/output_files)
DOT=$(wildcard $(MODELREP)*/output_files/*.dot)
PDF=$(DOT:%.dot=%.pdf)

DOTENGINE=dot
OPTION?=matlab

.PRECIOUS: $(PDF) $(DOT) $(OUTPUTREP)

all: clean $(EXAMPLES)
	@make -C . pdf 

models: $(EXAMPLES)
pdf: $(PDF)

%/output_files:
	@mkdir $@

%: $(MODELREP)% $(MODELREP)%/output_files
	@rm -rf $(MAINREP)$@_tmp.m 
	@echo "frontend('$@')" > $(MAINREP)$@_tmp.m  
	cd $(MAINREP) ; $(OPTION) < $@_tmp.m 
	@rm $(MAINREP)$@_tmp.m

%.pdf: %.dot
	$(DOTENGINE) -Tpdf $< -o $@

clean:
	rm -rf $(OUTPUTREP)

help:
	@echo make: compile all the models and generate the corresponding pdf
	@echo make models: compile all models
	@echo make pdf: convert each dot file into a pdf file 
	@echo make prozone: compile the prozone model 
	@echo make unary_vs_binary: compile the model with several time scales
