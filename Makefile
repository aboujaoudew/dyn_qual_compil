MAINREP=main/
MODELREP=case_studies/
EXAMPLES=$(notdir $(wildcard $(MODELREP)*))
OUTPUTREP=$(wildcard $(MODELREP)*/output_files)
DOT=$(wildcard $(MODELREP)*/output_files/*.dot)
PDF=$(DOT:%.dot=%.pdf)

DOTENGINE=dot
MATLAB=matlab

.PRECIOUS: $(PDF) $(DOT) $(OUTPUTREP)

all: $(EXAMPLES)
	@make -C . pdf 

models: $(EXAMPLES)
pdf: $(PDF)

%/output_files:
	@mkdir $@

%: $(MODELREP)% $(MODELREP)%/output_files
	@cd $(MAINREP) ; rm -rf tmp.m ; echo "frontend('$@')" > tmp.m ; $(MATLAB) < tmp.m 

%.pdf: %.dot
	@$(DOTENGINE) -Tpdf $< -o $@

clean:
	@rm -rf $(OUTPUTREP)

help:
	@echo make: compile all the models and generate the corresponfing pdf
	@echo make models: compile all models
	@echo make pdf: convert each dot file into a pdf file 
	@echo make prozone: compile the prozone model 
	@echo make unary_vs_binary: compile the model with several time scales
