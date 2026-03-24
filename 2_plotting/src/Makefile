F = gfortran

opt = -O3

#-----
all:  ../bin/fgenpdtdis ../bin/fgenstrdiprup ../bin/fgenmrf ../bin/fgensnap ../bin/fgentpdt ../bin/fgenmpinfo ../bin/GetASCIIwave ../bin/fgenazidel2xy ../bin/fgensnap2 ../bin/fgensnap_y ../bin/reverse_mrf_v2 ../bin/fgenpdtdis_mf ../bin/fgensnap_y_mf ../bin/fgensnap_mf ../bin/fgenstrdiprup_mf

#-----
../bin/fgenpdtdis : gen_pdtdis.f90 sub.gen.f90 FPSPACK.FOR
		$(F) -o $(@) $^ $(opt)

../bin/fgenstrdiprup : gen_strdiprup.f90 sub.gen.f90 FPSPACK.FOR
		$(F) -o $(@) $^ $(opt)

../bin/fgenmrf : gen_mrf.f90 sub.gen.f90 FPSPACK.FOR
		$(F) -o $(@) $^ $(opt)

../bin/fgensnap : gen_snap.f90 sub.gen.f90 FPSPACK.FOR
		$(F) -o $(@) $^ $(opt)

../bin/fgentpdt : gen_tpdt.f90 sub.gen.f90 FPSPACK.FOR
		$(F) -o $(@) $^ $(opt)

../bin/fgenmpinfo : gen_mpinfo.f90 sub.gen.f90 FPSPACK.FOR
		$(F) -o $(@) $^ $(opt)

../bin/GetASCIIwave : GetASCIIwave.f90
		$(F) -o $(@) $^ $(opt)

../bin/fgenazidel2xy : gen_azidel2xy.f90
		$(F) -o $(@) $^ $(opt)

../bin/fgensnap2 : gen_snap2.f90 sub.gen.f90 FPSPACK.FOR
		$(F) -o $(@) $^ $(opt)

../bin/fgensnap_y : gen_snap_y.f90 sub.gen.f90 FPSPACK.FOR
		$(F) -o $(@) $^ $(opt)

../bin/reverse_mrf_v2 : reverse_mrf_cli.f90 sub.gen.f90 FPSPACK.FOR
		$(F) -o $(@) $^ $(opt)

../bin/fgenpdtdis_mf : gen_pdtdis_mf.f90 sub.gen.f90 FPSPACK.FOR
		$(F) -o $(@) $^ $(opt)

../bin/fgensnap_y_mf : gen_snap_ymf.f90 sub.gen.f90 FPSPACK.FOR
		$(F) -o $(@) $^ $(opt)

../bin/fgensnap_mf : gen_snap_mf.f90 sub.gen.f90 FPSPACK.FOR
		$(F) -o $(@) $^ $(opt)

../bin/fgenstrdiprup_mf : gen_strdiprup_mf.f90 sub.gen.f90 FPSPACK.FOR
		$(F) -o $(@) $^ $(opt)

clean:
		rm ../bin/fgenpdtdis ../bin/fgenstrdiprup ../bin/fgenmrf ../bin/fgensnap ../bin/fgentpdt ../bin/fgenmpinfo ../bin/GetASCIIwave ../bin/fgenazidel2xy ../bin/fgensnap2 ../bin/fgensnap_y ../bin/reverse_mrf_v2 ../bin/fgenpdtdis_mf ../bin/fgensnap_y_mf ../bin/fgensnap_mf ../bin/fgenstrdiprup_mf
