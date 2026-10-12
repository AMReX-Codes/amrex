// Body of parser_exe_eval and parser_exe_eval_threaded.  Included inside
// those functions with `p` (bytecode) and `x` (variables) in scope.
// AMREX_PARSER_THREADED_BODY selects switch (0) or threaded (1) dispatch.

    // The top of the evaluation stack is kept in `top`.  A push first
    // spills `top` into pdata[psize++], so stack entry k lives in
    // pdata[k+1] (pdata[0] holds an unused value) and psize equals the
    // logical stack depth.  Local variables are stack entries.  They are
    // always spilled before they are read, because every statement
    // following an assignment starts with a push.
    double pdata[AMREX_PARSER_STACK_SIZE];
    int psize = 0;
    double top = 0.0;

#if AMREX_PARSER_THREADED_BODY
    // Threaded dispatch: each instruction jumps directly to the next one.
    static void* const dispatch_table[] = {
#define AMREX_PARSER_EXE_LABEL(name) &&L_PARSER_EXE_##name,
        AMREX_PARSER_EXE_TYPES(AMREX_PARSER_EXE_LABEL)
#undef AMREX_PARSER_EXE_LABEL
    };
#define AMREX_PARSER_CASE(X) L_##X:
#define AMREX_PARSER_NEXT goto *dispatch_table[*((parser_exe_t const*)p)]
    AMREX_PARSER_NEXT;
    {
        AMREX_PARSER_CASE(PARSER_EXE_NULL)
        {
            return top;
        }
#else
#define AMREX_PARSER_CASE(X) case X:
#define AMREX_PARSER_NEXT break
    while (*((parser_exe_t*)p) != PARSER_EXE_NULL) { // NOLINT
        switch (*((parser_exe_t*)p))
        {
#endif
        AMREX_PARSER_CASE(PARSER_EXE_NUMBER)
        {
            pdata[psize++] = top; // spill before reading locals
            top = ((ParserExeNumber*)p)->v;
            p += sizeof(ParserExeNumber);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_SYMBOL)
        {
            pdata[psize++] = top; // spill before reading locals
            int i = ((ParserExeSymbol*)p)->i;
            top = AMREX_PARSER_GET_DATA(i);
            p += sizeof(ParserExeSymbol);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_ADD)
        {
            double a = pdata[--psize]; // NOLINT
            top = a + top;
            p += sizeof(ParserExeADD);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_SUB_F)
        {
            double a = pdata[--psize]; // NOLINT
            top = a - top;
            p += sizeof(ParserExeSUB_F);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_SUB_B)
        {
            double a = pdata[--psize]; // NOLINT
            top = top - a;
            p += sizeof(ParserExeSUB_B);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_MUL)
        {
            double a = pdata[--psize]; // NOLINT
            top = a * top;
            p += sizeof(ParserExeMUL);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_DIV_F)
        {
            double a = pdata[--psize]; // NOLINT
            top = a / top;
            p += sizeof(ParserExeDIV_F);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_DIV_B)
        {
            double a = pdata[--psize]; // NOLINT
            top = top / a;
            p += sizeof(ParserExeDIV_B);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_F1)
        {
            top = parser_call_f1(((ParserExeF1*)p)->ftype, top);
            p += sizeof(ParserExeF1);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_F2_F)
        {
            double a = pdata[--psize]; // NOLINT
            top = parser_call_f2(((ParserExeF2_F*)p)->ftype, a, top);
            p += sizeof(ParserExeF2_F);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_F2_B)
        {
            double a = pdata[--psize]; // NOLINT
            top = parser_call_f2(((ParserExeF2_B*)p)->ftype, top, a);
            p += sizeof(ParserExeF2_B);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_ADD_VP)
        {
            pdata[psize++] = top; // spill before reading locals
            int i = ((ParserExeADD_VP*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            top = ((ParserExeADD_VP*)p)->v + d;
            p += sizeof(ParserExeADD_VP);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_SUB_VP)
        {
            pdata[psize++] = top; // spill before reading locals
            int i = ((ParserExeSUB_VP*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            top = ((ParserExeSUB_VP*)p)->v - d;
            p += sizeof(ParserExeSUB_VP);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_MUL_VP)
        {
            pdata[psize++] = top; // spill before reading locals
            int i = ((ParserExeMUL_VP*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            top = ((ParserExeMUL_VP*)p)->v * d;
            p += sizeof(ParserExeMUL_VP);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_DIV_VP)
        {
            pdata[psize++] = top; // spill before reading locals
            int i = ((ParserExeDIV_VP*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            top = ((ParserExeDIV_VP*)p)->v / d;
            p += sizeof(ParserExeDIV_VP);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_ADD_PP)
        {
            pdata[psize++] = top; // spill before reading locals
            int i = ((ParserExeADD_PP*)p)->i1;
            double d1 = AMREX_PARSER_GET_DATA(i);
            i     = ((ParserExeADD_PP*)p)->i2;
            double d2 = AMREX_PARSER_GET_DATA(i);
            top = d1+d2;
            p += sizeof(ParserExeADD_PP);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_SUB_PP)
        {
            pdata[psize++] = top; // spill before reading locals
            int i = ((ParserExeSUB_PP*)p)->i1;
            double d1 = AMREX_PARSER_GET_DATA(i);
            i     = ((ParserExeSUB_PP*)p)->i2;
            double d2 = AMREX_PARSER_GET_DATA(i);
            top = d1-d2;
            p += sizeof(ParserExeSUB_PP);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_MUL_PP)
        {
            pdata[psize++] = top; // spill before reading locals
            int i = ((ParserExeMUL_PP*)p)->i1;
            double d1 = AMREX_PARSER_GET_DATA(i);
            i     = ((ParserExeMUL_PP*)p)->i2;
            double d2 = AMREX_PARSER_GET_DATA(i);
            top = d1*d2;
            p += sizeof(ParserExeMUL_PP);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_DIV_PP)
        {
            pdata[psize++] = top; // spill before reading locals
            int i = ((ParserExeDIV_PP*)p)->i1;
            double d1 = AMREX_PARSER_GET_DATA(i);
            i     = ((ParserExeDIV_PP*)p)->i2;
            double d2 = AMREX_PARSER_GET_DATA(i);
            top = d1/d2;
            p += sizeof(ParserExeDIV_PP);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_ADD_VN)
        {
            top += ((ParserExeADD_VN*)p)->v;
            p += sizeof(ParserExeADD_VN);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_SUB_VN)
        {
            top = ((ParserExeSUB_VN*)p)->v - top;
            p += sizeof(ParserExeSUB_VN);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_MUL_VN)
        {
            top *= ((ParserExeMUL_VN*)p)->v;
            p += sizeof(ParserExeMUL_VN);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_DIV_VN)
        {
            top = ((ParserExeDIV_VN*)p)->v / top;
            p += sizeof(ParserExeDIV_VN);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_ADD_PN)
        {
            int i = ((ParserExeADD_PN*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            top += d;
            p += sizeof(ParserExeADD_PN);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_SUB_PN)
        {
            int i = ((ParserExeSUB_PN*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            if (((ParserExeSUB_PN*)p)->reverse) {
                top = top - d;
            } else {
                top = d - top;
            }
            p += sizeof(ParserExeSUB_PN);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_MUL_PN)
        {
            int i = ((ParserExeMUL_PN*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            top *= d;
            p += sizeof(ParserExeMUL_PN);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_DIV_PN)
        {
            int i = ((ParserExeDIV_PN*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            if (((ParserExeDIV_PN*)p)->reverse) {
                top = top / d;
            } else {
                top = d / top;
            }
            p += sizeof(ParserExeDIV_PN);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_SQUARE)
        {
            top *= top;
            p += sizeof(ParserExeSquare);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_POWI)
        {
            top = parser_math_powi(top, ((ParserExePOWI*)p)->i);
            p += sizeof(ParserExePOWI);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_FMA_VPV)
        {
            pdata[psize++] = top; // spill before reading locals
            int i = ((ParserExeFMA_VPV*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            top = ((ParserExeFMA_VPV*)p)->a * d + ((ParserExeFMA_VPV*)p)->b;
            p += sizeof(ParserExeFMA_VPV);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_IF)
        {
            double cond = top;
            top = pdata[--psize]; // NOLINT
            if (cond == 0.0) { // false branch
                p += ((ParserExeIF*)p)->offset;
            }
            p += sizeof(ParserExeIF);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_JUMP)
        {
            int offset = ((ParserExeJUMP*)p)->offset;
            p += sizeof(ParserExeJUMP) + offset;
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_USER_FN)
        {
            psize -= ((ParserExeUserFn*)p)->argc - 1;
            top = parser_detail::call_user_fn((ParserExeUserFn*)p, pdata+psize, top);
            p += sizeof(ParserExeUserFn);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_F1_P)
        {
            pdata[psize++] = top; // spill before reading locals
            int i = ((ParserExeF1_P*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            top = parser_call_f1(((ParserExeF1_P*)p)->ftype, d);
            p += sizeof(ParserExeF1_P);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_SQUARE_P)
        {
            pdata[psize++] = top; // spill before reading locals
            int i = ((ParserExeSquare_P*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            top = d*d;
            p += sizeof(ParserExeSquare_P);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_POWI_P)
        {
            pdata[psize++] = top; // spill before reading locals
            int i = ((ParserExePOWI_P*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            top = parser_math_powi(d, ((ParserExePOWI_P*)p)->n);
            p += sizeof(ParserExePOWI_P);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_F2_PP)
        {
            pdata[psize++] = top; // spill before reading locals
            int i = ((ParserExeF2_PP*)p)->i1;
            double d1 = AMREX_PARSER_GET_DATA(i);
            i     = ((ParserExeF2_PP*)p)->i2;
            double d2 = AMREX_PARSER_GET_DATA(i);
            top = parser_call_f2(((ParserExeF2_PP*)p)->ftype, d1, d2);
            p += sizeof(ParserExeF2_PP);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_F2_PV)
        {
            pdata[psize++] = top; // spill before reading locals
            int i = ((ParserExeF2_PV*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            top = parser_call_f2(((ParserExeF2_PV*)p)->ftype, d, ((ParserExeF2_PV*)p)->v);
            p += sizeof(ParserExeF2_PV);
            AMREX_PARSER_NEXT;
        }
        AMREX_PARSER_CASE(PARSER_EXE_F2_VP)
        {
            pdata[psize++] = top; // spill before reading locals
            int i = ((ParserExeF2_VP*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            top = parser_call_f2(((ParserExeF2_VP*)p)->ftype, ((ParserExeF2_VP*)p)->v, d);
            p += sizeof(ParserExeF2_VP);
            AMREX_PARSER_NEXT;
        }
#if AMREX_PARSER_THREADED_BODY
    }
#else
        default:
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,"parser_exe_eval: unknown node type");
        }
    }
    return top;
#endif
