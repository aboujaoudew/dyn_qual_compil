MAINREP=main/
MODELREP=case_studies/
EXAMPLES=$(notdir $(wildcard $(MODELREP)*))
OUTPUTREP=$(wildcard $(MODELREP)*/output_files)
DOT=$(wildcard $(MODELREP)*/output_files/*.dot)
PDF=$(DOT:%.dot=%.pdf)

DOTENGINE=dot
OPTION?=matlab
PROZONE=$(notdir $(wildcard $(MODELREP)*prozone*))
UNARY=$(notdir $(wildcard $(MODELREP)*unary_vs_binary*))

.PRECIOUS: $(PDF) $(DOT) $(OUTPUTREP)

all: clean $(EXAMPLES)
	@make -C . pdf 

models: $(EXAMPLES)
prozone: $(PROZONE)
unary_vs_binary: $(UNARY)

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
	@echo make Fig2_prozone242000: to generate the results for the prozone case study, with initial state 242000
	@echo make Fig3_prozone212000: to generate the results for the prozone case study, with initial state 212000
	@echo make Fig5a_unary_vs_binary200: to generate the results for the unary vs binary case study with initial state 200
	@echo make Fig5b_unary_vs_binary700: to generate the results for the unary vs binary case study with initial state 700
	@echo make prozone: compile the prozone model 
	@echo make unary_vs_binary: compile the model with several time scales
