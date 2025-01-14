using CxxWrap

const libpath = "build/libjuliaTPSA.so"
@wrapmodule(libpath)

da_init(3, 2, 400, false)

da = base(1)
x = 1 + da
x.print()