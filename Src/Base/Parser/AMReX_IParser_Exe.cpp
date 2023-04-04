#include <AMReX_IParser_Exe.H>

namespace amrex {

static int iparser_local_symbol_index (struct iparser_symbol* sym, Vector<char*>& local_variables)
{
    auto r = std::find_if(local_variables.rbegin(), local_variables.rend(),
                          [=] (char* i) { return std::strcmp(sym->name, i) == 0; });
    if (r != local_variables.rend()) {
        return static_cast<int>(std::distance(r, local_variables.rend())) - 1;
    } else {
        return -1;
    }
}

void
iparser_compile_exe_size (struct iparser_node* node, char*& p, std::size_t& exe_size,
                          int& max_stack_size, int& stack_size, Vector<char*>& local_variables)
{
    // In iparser_exe_eval, we push to the stack for NUMBER, SYMBOL, VP, PV, PP, and NEG_P.
    // In iparser_exe_eval, we pop the stack for ADD, SUB, MUL, DIV, F2, and IF.

    switch (node->type)
    {
    case IPARSER_NUMBER:
    {
        if (p) {
            auto *t = new(p) IParserExeNumber;
            p     += sizeof(IParserExeNumber);
            t->v = ((struct iparser_number*)node)->value;
        }
        exe_size += sizeof(IParserExeNumber);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case IPARSER_SYMBOL:
    {
        if (p) {
            auto *t = new(p) IParserExeSymbol;
            p     += sizeof(IParserExeSymbol);
            int lidx = iparser_local_symbol_index((struct iparser_symbol*)node, local_variables);
            if (lidx >= 0) {
                t->i = AMREX_IPARSER_LOCAL_IDX0 + lidx;
            } else {
                t->i = ((struct iparser_symbol*)node)->ip;
                if (t->i < 0) {
                    throw std::runtime_error(std::string("Unknown variable ")
                                             + ((struct iparser_symbol*)node)->name);
                }
            }
        }
        exe_size += sizeof(IParserExeSymbol);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case IPARSER_ADD:
    {
        if (node->l->type == IPARSER_NUMBER)
        {
            iparser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                     local_variables);
            if (p) {
                auto *t = new(p) IParserExeADD_VN;
                p     += sizeof(IParserExeADD_VN);
                t->v = ((struct iparser_number*)(node->l))->value;
            }
            exe_size += sizeof(IParserExeADD_VN);
        }
        else if (node->r->type == IPARSER_NUMBER)
        {
            iparser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                     local_variables);
            if (p) {
                auto *t = new(p) IParserExeADD_VN;
                p     += sizeof(IParserExeADD_VN);
                t->v = ((struct iparser_number*)(node->r))->value;
            }
            exe_size += sizeof(IParserExeADD_VN);
        }
        else if (node->l->type == IPARSER_SYMBOL)
        {
            iparser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                     local_variables);
            if (p) {
                auto *t = new(p) IParserExeADD_PN;
                p     += sizeof(IParserExeADD_PN);
                int lidx = iparser_local_symbol_index((struct iparser_symbol*)(node->l),
                                                      local_variables);
                if (lidx >= 0) {
                    t->i = AMREX_IPARSER_LOCAL_IDX0 + lidx;
                } else {
                    t->i = ((struct iparser_symbol*)(node->l))->ip;
                    if (t->i < 0) {
                        throw std::runtime_error(std::string("Unknown variable ")
                                                 + ((struct iparser_symbol*)node->l)->name);
                    }
                }
            }
            exe_size += sizeof(IParserExeADD_PN);
        }
        else if (node->r->type == IPARSER_SYMBOL)
        {
            iparser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                     local_variables);
            if (p) {
                auto *t = new(p) IParserExeADD_PN;
                p     += sizeof(IParserExeADD_PN);
                int lidx = iparser_local_symbol_index((struct iparser_symbol*)(node->r),
                                                      local_variables);
                if (lidx >= 0) {
                    t->i = AMREX_IPARSER_LOCAL_IDX0 + lidx;
                } else {
                    t->i = ((struct iparser_symbol*)(node->r))->ip;
                    if (t->i < 0) {
                        throw std::runtime_error(std::string("Unknown variable ")
                                                 + ((struct iparser_symbol*)node->r)->name);
                    }
                }
            }
            exe_size += sizeof(IParserExeADD_PN);
        }
        else
        {
            int d1 = iparser_ast_depth(node->l);
            int d2 = iparser_ast_depth(node->r);
            if (d1 < d2) {
                iparser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                         local_variables);
                iparser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                         local_variables);
            } else {
                iparser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                         local_variables);
                iparser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                         local_variables);
            }
            if (p) {
                new(p)      IParserExeADD;
                p += sizeof(IParserExeADD);
            }
            exe_size += sizeof(IParserExeADD);
            --stack_size;
        }
        // no need to update max_stack_size
        break;
    }
    case IPARSER_SUB:
    {
        if (node->l->type == IPARSER_NUMBER)
        {
            iparser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                     local_variables);
            if (p) {
                auto *t = new(p) IParserExeSUB_VN;
                p     += sizeof(IParserExeSUB_VN);
                t->v = ((struct iparser_number*)(node->l))->value;
            }
            exe_size += sizeof(IParserExeSUB_VN);
        }
        else if (node->r->type == IPARSER_NUMBER)
        {
            iparser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                     local_variables);
            if (p) {
                auto *t = new(p) IParserExeADD_VN;
                p     += sizeof(IParserExeADD_VN);
                t->v = -(((struct iparser_number*)(node->r))->value);
            }
            exe_size += sizeof(IParserExeADD_VN);
        }
        else if (node->l->type == IPARSER_SYMBOL)
        {
            iparser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                     local_variables);
            if (p) {
                auto *t = new(p) IParserExeSUB_PN;
                p     += sizeof(IParserExeSUB_PN);
                int lidx = iparser_local_symbol_index((struct iparser_symbol*)(node->l),
                                                      local_variables);
                if (lidx >= 0) {
                    t->i = AMREX_IPARSER_LOCAL_IDX0 + lidx;
                } else {
                    t->i = ((struct iparser_symbol*)(node->l))->ip;
                    if (t->i < 0) {
                        throw std::runtime_error(std::string("Unknown variable ")
                                                 + ((struct iparser_symbol*)node->l)->name);
                    }
                }
                t->sign = 1.0;
            }
            exe_size += sizeof(IParserExeSUB_PN);
        }
        else if (node->r->type == IPARSER_SYMBOL)
        {
            iparser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                     local_variables);
            if (p) {
                auto *t = new(p) IParserExeSUB_PN;
                p     += sizeof(IParserExeSUB_PN);
                int lidx = iparser_local_symbol_index((struct iparser_symbol*)(node->r),
                                                      local_variables);
                if (lidx >= 0) {
                    t->i = AMREX_IPARSER_LOCAL_IDX0 + lidx;
                } else {
                    t->i = ((struct iparser_symbol*)(node->r))->ip;
                    if (t->i < 0) {
                        throw std::runtime_error(std::string("Unknown variable ")
                                                 + ((struct iparser_symbol*)node->r)->name);
                    }
                }
                t->sign = -1.0;
            }
            exe_size += sizeof(IParserExeSUB_PN);
        }
        else
        {
            int d1 = iparser_ast_depth(node->l);
            int d2 = iparser_ast_depth(node->r);
            if (d1 < d2) {
                iparser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                         local_variables);
                iparser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                         local_variables);
            } else {
                iparser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                         local_variables);
                iparser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                         local_variables);
            }
            if (p) {
                auto *t = new(p) IParserExeSUB;
                p     += sizeof(IParserExeSUB);
                t->sign = (d1 < d2) ? -1.0 : 1.0;
            }
            exe_size += sizeof(IParserExeSUB);
            --stack_size;
        }
        // no need to update max_stack_size
        break;
    }
    case IPARSER_MUL:
    {
        if (node->l->type == IPARSER_NUMBER)
        {
            iparser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                     local_variables);
            if (p) {
                auto *t = new(p) IParserExeMUL_VN;
                p     += sizeof(IParserExeMUL_VN);
                t->v = ((struct iparser_number*)(node->l))->value;
            }
            exe_size += sizeof(IParserExeMUL_VN);
        }
        else if (node->r->type == IPARSER_NUMBER)
        {
            iparser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                     local_variables);
            if (p) {
                auto *t = new(p) IParserExeMUL_VN;
                p     += sizeof(IParserExeMUL_VN);
                t->v = ((struct iparser_number*)(node->r))->value;
            }
            exe_size += sizeof(IParserExeMUL_VN);
        }
        else if (node->l->type == IPARSER_SYMBOL)
        {
            iparser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                     local_variables);
            if (p) {
                auto *t = new(p) IParserExeMUL_PN;
                p     += sizeof(IParserExeMUL_PN);
                int lidx = iparser_local_symbol_index((struct iparser_symbol*)(node->l),
                                                      local_variables);
                if (lidx >= 0) {
                    t->i = AMREX_IPARSER_LOCAL_IDX0 + lidx;
                } else {
                    t->i = ((struct iparser_symbol*)(node->l))->ip;
                    if (t->i < 0) {
                        throw std::runtime_error(std::string("Unknown variable ")
                                                 + ((struct iparser_symbol*)node->l)->name);
                    }
                }
            }
            exe_size += sizeof(IParserExeMUL_PN);
        }
        else if (node->r->type == IPARSER_SYMBOL)
        {
            iparser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                     local_variables);
            if (p) {
                auto *t = new(p) IParserExeMUL_PN;
                p     += sizeof(IParserExeMUL_PN);
                int lidx = iparser_local_symbol_index((struct iparser_symbol*)(node->r),
                                                      local_variables);
                if (lidx >= 0) {
                    t->i = AMREX_IPARSER_LOCAL_IDX0 + lidx;
                } else {
                    t->i = ((struct iparser_symbol*)(node->r))->ip;
                    if (t->i < 0) {
                        throw std::runtime_error(std::string("Unknown variable ")
                                                 + ((struct iparser_symbol*)node->r)->name);
                    }
                }
            }
            exe_size += sizeof(IParserExeMUL_PN);
        }
        else
        {
            int d1 = iparser_ast_depth(node->l);
            int d2 = iparser_ast_depth(node->r);
            if (d1 < d2) {
                iparser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                         local_variables);
                iparser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                         local_variables);
            } else {
                iparser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                         local_variables);
                iparser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                         local_variables);
            }
            if (p) {
                new(p)      IParserExeMUL;
                p += sizeof(IParserExeMUL);
            }
            exe_size += sizeof(IParserExeMUL);
            --stack_size;
        }
        // no need to update max_stack_size
        break;
    }
    case IPARSER_DIV:
    {
        if (node->l->type == IPARSER_NUMBER)
        {
            iparser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                     local_variables);
            if (p) {
                auto *t = new(p) IParserExeDIV_VN;
                p     += sizeof(IParserExeDIV_VN);
                t->v = ((struct iparser_number*)(node->l))->value;
            }
            exe_size += sizeof(IParserExeDIV_VN);
        }
        else if (node->r->type == IPARSER_NUMBER)
        {
            iparser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                     local_variables);
            if (p) {
                auto *t = new(p) IParserExeDIV_NV;
                p     += sizeof(IParserExeDIV_NV);
                t->v = ((struct iparser_number*)(node->r))->value;
            }
            exe_size += sizeof(IParserExeDIV_NV);
        }
        else if (node->l->type == IPARSER_SYMBOL)
        {
            iparser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                    local_variables);
            if (p) {
                auto *t = new(p) IParserExeDIV_PN;
                p     += sizeof(IParserExeDIV_PN);
                int lidx = iparser_local_symbol_index((struct iparser_symbol*)(node->l),
                                                      local_variables);
                if (lidx >= 0) {
                    t->i = AMREX_IPARSER_LOCAL_IDX0 + lidx;
                } else {
                    t->i = ((struct iparser_symbol*)(node->l))->ip;
                    if (t->i < 0) {
                        throw std::runtime_error(std::string("Unknown variable ")
                                                 + ((struct iparser_symbol*)node->l)->name);
                    }
                }
                t->reverse = false;
            }
            exe_size += sizeof(IParserExeDIV_PN);
        }
        else if (node->r->type == IPARSER_SYMBOL)
        {
            iparser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                     local_variables);
            if (p) {
                auto *t = new(p) IParserExeDIV_PN;
                p     += sizeof(IParserExeDIV_PN);
                int lidx = iparser_local_symbol_index((struct iparser_symbol*)(node->r),
                                                      local_variables);
                if (lidx >= 0) {
                    t->i = AMREX_IPARSER_LOCAL_IDX0 + lidx;
                } else {
                    t->i = ((struct iparser_symbol*)(node->r))->ip;
                    if (t->i < 0) {
                        throw std::runtime_error(std::string("Unknown variable ")
                                                 + ((struct iparser_symbol*)node->r)->name);
                    }
                }
                t->reverse = true;
            }
            exe_size += sizeof(IParserExeDIV_PN);
        }
        else
        {
            int d1 = iparser_ast_depth(node->l);
            int d2 = iparser_ast_depth(node->r);
            if (d1 < d2) {
                iparser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                         local_variables);
                iparser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                         local_variables);
                if (p) {
                    new(p)      IParserExeDIV_B;
                    p += sizeof(IParserExeDIV_B);
                }
                exe_size += sizeof(IParserExeDIV_B);
            } else {
                iparser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                         local_variables);
                iparser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                         local_variables);
                if (p) {
                    new(p)      IParserExeDIV_F;
                    p += sizeof(IParserExeDIV_F);
                }
                exe_size += sizeof(IParserExeDIV_F);
            }
            --stack_size;
        }
        // no need to update max_grid_size
        break;
    }
    case IPARSER_NEG:
    {
        iparser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size, local_variables);
        if (p) {
            new(p)      IParserExeNEG;
            p += sizeof(IParserExeNEG);
        }
        exe_size += sizeof(IParserExeNEG);
        break;
    }
    case IPARSER_F1:
    {
        iparser_compile_exe_size(((struct iparser_f1*)node)->l,
                                 p, exe_size, max_stack_size, stack_size, local_variables);
        if (p) {
            auto *t = new(p) IParserExeF1;
            p     += sizeof(IParserExeF1);
            t->ftype = ((struct iparser_f1*)node)->ftype;
        }
        exe_size += sizeof(IParserExeF1);
        break;
    }
    case IPARSER_F2:
    {
        int d1 = iparser_ast_depth(((struct iparser_f2*)node)->l);
        int d2 = iparser_ast_depth(((struct iparser_f2*)node)->r);
        if (d1 < d2) {
            iparser_compile_exe_size(((struct iparser_f2*)node)->r,
                                     p, exe_size, max_stack_size, stack_size, local_variables);
            iparser_compile_exe_size(((struct iparser_f2*)node)->l,
                                     p, exe_size, max_stack_size, stack_size, local_variables);
            if (p) {
                auto *t = new(p) IParserExeF2_B;
                p     += sizeof(IParserExeF2_B);
                t->ftype = ((struct iparser_f2*)node)->ftype;
            }
            exe_size += sizeof(IParserExeF2_B);
        } else {
            iparser_compile_exe_size(((struct iparser_f2*)node)->l,
                                     p, exe_size, max_stack_size, stack_size, local_variables);
            iparser_compile_exe_size(((struct iparser_f2*)node)->r,
                                     p, exe_size, max_stack_size, stack_size, local_variables);
            if (p) {
                auto *t = new(p) IParserExeF2_F;
                p     += sizeof(IParserExeF2_F);
                t->ftype = ((struct iparser_f2*)node)->ftype;
            }
            exe_size += sizeof(IParserExeF2_F);
        }
        --stack_size;
        break;
    }
    case IPARSER_F3:
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(((struct iparser_f3*)node)->ftype == IPARSER_IF,
                                         "iparser_compile: unknown f3 type");
        iparser_compile_exe_size(((struct iparser_f3*)node)->n1,
                                 p, exe_size, max_stack_size, stack_size, local_variables);

        IParserExeIF* tif = nullptr;
        char* psave = nullptr;
        if (p) {
            tif = new(p) IParserExeIF;
            p  += sizeof(IParserExeIF);
            psave = p;
        }
        exe_size += sizeof(IParserExeIF);
        --stack_size;
        auto stack_size_save = stack_size;

        iparser_compile_exe_size(((struct iparser_f3*)node)->n2,
                                 p, exe_size, max_stack_size, stack_size, local_variables);

        IParserExeJUMP* tjump = nullptr;
        if (p) {
            tjump = new(p) IParserExeJUMP;
            p +=    sizeof(IParserExeJUMP);
        }
        exe_size += sizeof(IParserExeJUMP);

        if (psave) {
            tif->offset = static_cast<int>(p-psave);
        }
        stack_size = stack_size_save;

        psave = p;
        iparser_compile_exe_size(((struct iparser_f3*)node)->n3,
                                 p, exe_size, max_stack_size, stack_size, local_variables);
        if (tjump) {
            tjump->offset = static_cast<int>(p-psave);
        }

        break;
    }
    case IPARSER_ASSIGN:
    {
        auto *asgn = (struct iparser_assign*)node;
        local_variables.push_back(asgn->s->name);
        iparser_compile_exe_size(asgn->v, p, exe_size, max_stack_size, stack_size, local_variables);
        break;
    }
    case IPARSER_LIST:
    {
        iparser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size, local_variables);
        iparser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size, local_variables);
        break;
    }
    case IPARSER_ADD_VP:
    {
        if (p) {
            auto *t = new(p) IParserExeADD_VP;
            p     += sizeof(IParserExeADD_VP);
            int lidx = iparser_local_symbol_index((struct iparser_symbol*)(node->r), local_variables);
            if (lidx >= 0) {
                t->i = AMREX_IPARSER_LOCAL_IDX0 + lidx;
            } else {
                t->i = node->rip;
            }
            t->v = node->lvp.v;
        }
        exe_size += sizeof(IParserExeADD_VP);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case IPARSER_SUB_VP:
    {
        if (p) {
            auto *t = new(p) IParserExeSUB_VP;
            p     += sizeof(IParserExeSUB_VP);
            int lidx = iparser_local_symbol_index((struct iparser_symbol*)(node->r), local_variables);
            if (lidx >= 0) {
                t->i = AMREX_IPARSER_LOCAL_IDX0 + lidx;
            } else {
                t->i = node->rip;
            }
            t->v = node->lvp.v;
        }
        exe_size += sizeof(IParserExeSUB_VP);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case IPARSER_MUL_VP:
    {
        if (p) {
            auto *t = new(p) IParserExeMUL_VP;
            p     += sizeof(IParserExeMUL_VP);
            int lidx = iparser_local_symbol_index((struct iparser_symbol*)(node->r), local_variables);
            if (lidx >= 0) {
                t->i = AMREX_IPARSER_LOCAL_IDX0 + lidx;
            } else {
                t->i = node->rip;
            }
            t->v = node->lvp.v;
        }
        exe_size += sizeof(IParserExeMUL_VP);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case IPARSER_DIV_VP:
    {
        if (p) {
            auto *t = new(p) IParserExeDIV_VP;
            p     += sizeof(IParserExeDIV_VP);
            int lidx = iparser_local_symbol_index((struct iparser_symbol*)(node->r), local_variables);
            if (lidx >= 0) {
                t->i = AMREX_IPARSER_LOCAL_IDX0 + lidx;
            } else {
                t->i = node->rip;
            }
            t->v = node->lvp.v;
        }
        exe_size += sizeof(IParserExeDIV_VP);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case IPARSER_DIV_PV:
    {
        if (p) {
            auto *t = new(p) IParserExeDIV_PV;
            p     += sizeof(IParserExeDIV_PV);
            int lidx = iparser_local_symbol_index((struct iparser_symbol*)(node->r), local_variables);
            if (lidx >= 0) {
                t->i = AMREX_IPARSER_LOCAL_IDX0 + lidx;
            } else {
                t->i = node->rip;
            }
            t->v = node->lvp.v;
        }
        exe_size += sizeof(IParserExeDIV_PV);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case IPARSER_ADD_PP:
    {
        if (p) {
            auto *t = new(p) IParserExeADD_PP;
            p     += sizeof(IParserExeADD_PP);
            int li1 = iparser_local_symbol_index((struct iparser_symbol*)(node->l), local_variables);
            int li2 = iparser_local_symbol_index((struct iparser_symbol*)(node->r), local_variables);
            if (li1 >= 0) {
                t->i1 = AMREX_IPARSER_LOCAL_IDX0 + li1;
            } else {
                t->i1 = node->lvp.ip;
            }
            if (li2 >= 0) {
                t->i2 = AMREX_IPARSER_LOCAL_IDX0 + li2;
            } else {
                t->i2 = node->rip;
            }
        }
        exe_size += sizeof(IParserExeADD_PP);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case IPARSER_SUB_PP:
    {
        if (p) {
            auto *t = new(p) IParserExeSUB_PP;
            p     += sizeof(IParserExeSUB_PP);
            int li1 = iparser_local_symbol_index((struct iparser_symbol*)(node->l), local_variables);
            int li2 = iparser_local_symbol_index((struct iparser_symbol*)(node->r), local_variables);
            if (li1 >= 0) {
                t->i1 = AMREX_IPARSER_LOCAL_IDX0 + li1;
            } else {
                t->i1 = node->lvp.ip;
            }
            if (li2 >= 0) {
                t->i2 = AMREX_IPARSER_LOCAL_IDX0 + li2;
            } else {
                t->i2 = node->rip;
            }
        }
        exe_size += sizeof(IParserExeSUB_PP);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case IPARSER_MUL_PP:
    {
        if (p) {
            auto *t = new(p) IParserExeMUL_PP;
            p     += sizeof(IParserExeMUL_PP);
            int li1 = iparser_local_symbol_index((struct iparser_symbol*)(node->l), local_variables);
            int li2 = iparser_local_symbol_index((struct iparser_symbol*)(node->r), local_variables);
            if (li1 >= 0) {
                t->i1 = AMREX_IPARSER_LOCAL_IDX0 + li1;
            } else {
                t->i1 = node->lvp.ip;
            }
            if (li2 >= 0) {
                t->i2 = AMREX_IPARSER_LOCAL_IDX0 + li2;
            } else {
                t->i2 = node->rip;
            }
        }
        exe_size += sizeof(IParserExeMUL_PP);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case IPARSER_DIV_PP:
    {
        if (p) {
            auto *t = new(p) IParserExeDIV_PP;
            p     += sizeof(IParserExeDIV_PP);
            int li1 = iparser_local_symbol_index((struct iparser_symbol*)(node->l), local_variables);
            int li2 = iparser_local_symbol_index((struct iparser_symbol*)(node->r), local_variables);
            if (li1 >= 0) {
                t->i1 = AMREX_IPARSER_LOCAL_IDX0 + li1;
            } else {
                t->i1 = node->lvp.ip;
            }
            if (li2 >= 0) {
                t->i2 = AMREX_IPARSER_LOCAL_IDX0 + li2;
            } else {
                t->i2 = node->rip;
            }
        }
        exe_size += sizeof(IParserExeDIV_PP);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case IPARSER_NEG_P:
    {
        if (p) {
            auto *t = new(p) IParserExeNEG_P;
            p     += sizeof(IParserExeNEG_P);
            int lidx = iparser_local_symbol_index((struct iparser_symbol*)(node->l), local_variables);
            if (lidx >= 0) {
                t->i = AMREX_IPARSER_LOCAL_IDX0 + lidx;
            } else {
                t->i = node->lvp.ip;
            }
        }
        exe_size += sizeof(IParserExeNEG_P);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    default:
        amrex::Abort("iparser_compile: unknown node type " + std::to_string(node->type));
    }
}

}
