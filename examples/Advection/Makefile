help:
	@echo "make run"

run: example.mesh example.mtx example.sol.mtx

viz: example_viz.mesh example_viz.bb

#example.mesh:
#	../../bin/wmesh-build --nteta 13 --nradius 17 -o $@ --dim 2

example_viz.mesh:example.mesh
	../../bin/wmesh-refine $< -o $@ -d 2 -n lagrange -o $@ --dg
example_viz.bb:example.sol.mtx
	../../bin/wmesh-convert $< -o $@

%.mtx:%.mesh
	../../bin/wmesh-sparsedg $< -d 2 -n lagrange -o $*

%.sol.mtx:%.mtx
	../../../WLS/bin/wls-gesv $< $*_rhs.mtx  -o $@

clean:
	\rm -f *.mtx *.mesh *.bb *.sol



