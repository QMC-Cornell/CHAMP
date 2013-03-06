
BINDIR   = ../src/
INPUT    = pw2blip.in
INPUT0   = pw2blip.0.in
INPUT1   = pw2blip.1.in
INPUT2   = pw2blip.2.in
OUTPUT0  = pw2blip.0.out
OUTPUT1  = pw2blip.1.out
OUTPUT2  = pw2blip.2.out

default: 
	cp $(INPUT) $(INPUT0)
	$(BINDIR)/blip < $(INPUT0) > $(OUTPUT0)
	mv bwfn.data bwfn.data.0
	rm $(INPUT0)
	sed 's/ 0 / 1 /' $(INPUT) > $(INPUT1)
	$(BINDIR)/blip < $(INPUT1) > $(OUTPUT1)
	mv bwfn.data bwfn.data.1
	rm $(INPUT1)
	sed 's/ 0 / 2 /' $(INPUT) > $(INPUT2)
	$(BINDIR)/blip < $(INPUT2) > $(OUTPUT2)
	mv bwfn.data bwfn.data.2
	rm $(INPUT2)

clean:
	rm -rf bwfn.data* *.out
