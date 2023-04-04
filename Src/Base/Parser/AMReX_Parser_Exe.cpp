#include <AMReX_Parser_Exe.H>

namespace amrex {

static int parser_local_symbol_index (struct parser_symbol* sym, Vector<char*>& local_variables)
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
parser_compile_exe_size (struct parser_node* node, char*& p, std::size_t& exe_size,
                         int& max_stack_size, int& stack_size, Vector<char*>& local_variables)
{
    // In parser_exe_eval, we push to the stack for NUMBER, SYMBOL, VP, PP, and NEG_P.
    // In parser_exe_eval, we pop the stack for ADD, SUB, MUL, DIV, F2, and IF.

    switch (node->type)
    {
    case PARSER_NUMBER:
    {
        if (p) {
            auto *t = new(p) ParserExeNumber;
            p     += sizeof(ParserExeNumber);
            t->v = ((struct parser_number*)node)->value;
        }
        exe_size += sizeof(ParserExeNumber);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case PARSER_SYMBOL:
    {
        if (p) {
            auto *t = new(p) ParserExeSymbol;
            p     += sizeof(ParserExeSymbol);
            int lidx = parser_local_symbol_index((struct parser_symbol*)node, local_variables);
            if (lidx >= 0) {
                t->i = AMREX_PARSER_LOCAL_IDX0 + lidx;
            } else {
                t->i = ((struct parser_symbol*)node)->ip;
                if (t->i < 0) {
                    throw std::runtime_error(std::string("Unknown variable ")
                                             + ((struct parser_symbol*)node)->name);
                }
            }
        }
        exe_size += sizeof(ParserExeSymbol);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case PARSER_ADD:
    {
        if (node->l->type == PARSER_NUMBER)
        {
            parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                    local_variables);
            if (p) {
                auto *t = new(p) ParserExeADD_VN;
                p     += sizeof(ParserExeADD_VN);
                t->v = ((struct parser_number*)(node->l))->value;
            }
            exe_size += sizeof(ParserExeADD_VN);
        }
        else if (node->r->type == PARSER_NUMBER)
        {
            parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                    local_variables);
            if (p) {
                auto *t = new(p) ParserExeADD_VN;
                p     += sizeof(ParserExeADD_VN);
                t->v = ((struct parser_number*)(node->r))->value;
            }
            exe_size += sizeof(ParserExeADD_VN);
        }
        else if (node->l->type == PARSER_SYMBOL)
        {
            parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                    local_variables);
            if (p) {
                auto *t = new(p) ParserExeADD_PN;
                p     += sizeof(ParserExeADD_PN);
                int lidx = parser_local_symbol_index((struct parser_symbol*)(node->l),
                                                     local_variables);
                if (lidx >= 0) {
                    t->i = AMREX_PARSER_LOCAL_IDX0 + lidx;
                } else {
                    t->i = ((struct parser_symbol*)(node->l))->ip;
                    if (t->i < 0) {
                        throw std::runtime_error(std::string("Unknown variable ")
                                                 + ((struct parser_symbol*)node->l)->name);
                    }
                }
            }
            exe_size += sizeof(ParserExeADD_PN);
        }
        else if (node->r->type == PARSER_SYMBOL)
        {
            parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                    local_variables);
            if (p) {
                auto *t = new(p) ParserExeADD_PN;
                p     += sizeof(ParserExeADD_PN);
                int lidx = parser_local_symbol_index((struct parser_symbol*)(node->r),
                                                     local_variables);
                if (lidx >= 0) {
                    t->i = AMREX_PARSER_LOCAL_IDX0 + lidx;
                } else {
                    t->i = ((struct parser_symbol*)(node->r))->ip;
                    if (t->i < 0) {
                        throw std::runtime_error(std::string("Unknown variable ")
                                                 + ((struct parser_symbol*)node->r)->name);
                    }
                }
            }
            exe_size += sizeof(ParserExeADD_PN);
        }
        else
        {
            int d1 = parser_ast_depth(node->l);
            int d2 = parser_ast_depth(node->r);
            if (d1 < d2) {
                parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                        local_variables);
                parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                        local_variables);
            } else {
                parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                        local_variables);
                parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                        local_variables);
            }
            if (p) {
                new(p)      ParserExeADD;
                p += sizeof(ParserExeADD);
            }
            exe_size += sizeof(ParserExeADD);
            --stack_size;
        }
        // no need to update max_stack_size
        break;
    }
    case PARSER_SUB:
    {
        if (node->l->type == PARSER_NUMBER)
        {
            parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                    local_variables);
            if (p) {
                auto *t = new(p) ParserExeSUB_VN;
                p     += sizeof(ParserExeSUB_VN);
                t->v = ((struct parser_number*)(node->l))->value;
            }
            exe_size += sizeof(ParserExeSUB_VN);
        }
        else if (node->r->type == PARSER_NUMBER)
        {
            parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                    local_variables);
            if (p) {
                auto *t = new(p) ParserExeADD_VN;
                p     += sizeof(ParserExeADD_VN);
                t->v = -(((struct parser_number*)(node->r))->value);
            }
            exe_size += sizeof(ParserExeADD_VN);
        }
        else if (node->l->type == PARSER_SYMBOL)
        {
            parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                    local_variables);
            if (p) {
                auto *t = new(p) ParserExeSUB_PN;
                p     += sizeof(ParserExeSUB_PN);
                int lidx = parser_local_symbol_index((struct parser_symbol*)(node->l),
                                                     local_variables);
                if (lidx >= 0) {
                    t->i = AMREX_PARSER_LOCAL_IDX0 + lidx;
                } else {
                    t->i = ((struct parser_symbol*)(node->l))->ip;
                    if (t->i < 0) {
                        throw std::runtime_error(std::string("Unknown variable ")
                                                 + ((struct parser_symbol*)node->l)->name);
                    }
                }
                t->sign = 1.0;
            }
            exe_size += sizeof(ParserExeSUB_PN);
        }
        else if (node->r->type == PARSER_SYMBOL)
        {
            parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                    local_variables);
            if (p) {
                auto *t = new(p) ParserExeSUB_PN;
                p     += sizeof(ParserExeSUB_PN);
                int lidx = parser_local_symbol_index((struct parser_symbol*)(node->r),
                                                     local_variables);
                if (lidx >= 0) {
                    t->i = AMREX_PARSER_LOCAL_IDX0 + lidx;
                } else {
                    t->i = ((struct parser_symbol*)(node->r))->ip;
                    if (t->i < 0) {
                        throw std::runtime_error(std::string("Unknown variable ")
                                                 + ((struct parser_symbol*)node->r)->name);
                    }
                }
                t->sign = -1.0;
            }
            exe_size += sizeof(ParserExeSUB_PN);
        }
        else
        {
            int d1 = parser_ast_depth(node->l);
            int d2 = parser_ast_depth(node->r);
            if (d1 < d2) {
                parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                        local_variables);
                parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                        local_variables);
            } else {
                parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                        local_variables);
                parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                        local_variables);
            }
            if (p) {
                auto *t = new(p) ParserExeSUB;
                p     += sizeof(ParserExeSUB);
                t->sign = (d1 < d2) ? -1.0 : 1.0;
            }
            exe_size += sizeof(ParserExeSUB);
            --stack_size;
        }
        // no need to update max_stack_size
        break;
    }
    case PARSER_MUL:
    {
        if (node->l->type == PARSER_NUMBER)
        {
            parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                    local_variables);
            if (p) {
                auto *t = new(p) ParserExeMUL_VN;
                p     += sizeof(ParserExeMUL_VN);
                t->v = ((struct parser_number*)(node->l))->value;
            }
            exe_size += sizeof(ParserExeMUL_VN);
        }
        else if (node->r->type == PARSER_NUMBER)
        {
            parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                    local_variables);
            if (p) {
                auto *t = new(p) ParserExeMUL_VN;
                p     += sizeof(ParserExeMUL_VN);
                t->v = ((struct parser_number*)(node->r))->value;
            }
            exe_size += sizeof(ParserExeMUL_VN);
        }
        else if (node->l->type == PARSER_SYMBOL)
        {
            parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                    local_variables);
            if (p) {
                auto *t = new(p) ParserExeMUL_PN;
                p     += sizeof(ParserExeMUL_PN);
                int lidx = parser_local_symbol_index((struct parser_symbol*)(node->l),
                                                     local_variables);
                if (lidx >= 0) {
                    t->i = AMREX_PARSER_LOCAL_IDX0 + lidx;
                } else {
                    t->i = ((struct parser_symbol*)(node->l))->ip;
                    if (t->i < 0) {
                        throw std::runtime_error(std::string("Unknown variable ")
                                                 + ((struct parser_symbol*)node->l)->name);
                    }
                }
            }
            exe_size += sizeof(ParserExeMUL_PN);
        }
        else if (node->r->type == PARSER_SYMBOL)
        {
            parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                    local_variables);
            if (p) {
                auto *t = new(p) ParserExeMUL_PN;
                p     += sizeof(ParserExeMUL_PN);
                int lidx = parser_local_symbol_index((struct parser_symbol*)(node->r),
                                                     local_variables);
                if (lidx >= 0) {
                    t->i = AMREX_PARSER_LOCAL_IDX0 + lidx;
                } else {
                    t->i = ((struct parser_symbol*)(node->r))->ip;
                    if (t->i < 0) {
                        throw std::runtime_error(std::string("Unknown variable ")
                                                 + ((struct parser_symbol*)node->r)->name);
                    }
                }
            }
            exe_size += sizeof(ParserExeMUL_PN);
        }
        else
        {
            int d1 = parser_ast_depth(node->l);
            int d2 = parser_ast_depth(node->r);
            if (d1 < d2) {
                parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                        local_variables);
                parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                        local_variables);
            } else {
                parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                        local_variables);
                parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                        local_variables);
            }
            if (p) {
                new(p)      ParserExeMUL;
                p += sizeof(ParserExeMUL);
            }
            exe_size += sizeof(ParserExeMUL);
            --stack_size;
        }
        // no need to update max_stack_size
        break;
    }
    case PARSER_DIV:
    {
        if (node->l->type == PARSER_NUMBER)
        {
            parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                    local_variables);
            if (p) {
                auto *t = new(p) ParserExeDIV_VN;
                p     += sizeof(ParserExeDIV_VN);
                t->v = ((struct parser_number*)(node->l))->value;
            }
            exe_size += sizeof(ParserExeDIV_VN);
        }
        else if (node->r->type == PARSER_NUMBER)
        {
            parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                    local_variables);
            if (p) {
                auto *t = new(p) ParserExeMUL_VN;
                p     += sizeof(ParserExeMUL_VN);
                t->v = 1.0 / ((struct parser_number*)(node->r))->value;
            }
            exe_size += sizeof(ParserExeMUL_VN);
        }
        else if (node->l->type == PARSER_SYMBOL)
        {
            parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                    local_variables);
            if (p) {
                auto *t = new(p) ParserExeDIV_PN;
                p     += sizeof(ParserExeDIV_PN);
                int lidx = parser_local_symbol_index((struct parser_symbol*)(node->l),
                                                     local_variables);
                if (lidx >= 0) {
                    t->i = AMREX_PARSER_LOCAL_IDX0 + lidx;
                } else {
                    t->i = ((struct parser_symbol*)(node->l))->ip;
                    if (t->i < 0) {
                        throw std::runtime_error(std::string("Unknown variable ")
                                                 + ((struct parser_symbol*)node->l)->name);
                    }
                }
                t->reverse = false;
            }
            exe_size += sizeof(ParserExeDIV_PN);
        }
        else if (node->r->type == PARSER_SYMBOL)
        {
            parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                    local_variables);
            if (p) {
                auto *t = new(p) ParserExeDIV_PN;
                p     += sizeof(ParserExeDIV_PN);
                int lidx = parser_local_symbol_index((struct parser_symbol*)(node->r),
                                                     local_variables);
                if (lidx >= 0) {
                    t->i = AMREX_PARSER_LOCAL_IDX0 + lidx;
                } else {
                    t->i = ((struct parser_symbol*)(node->r))->ip;
                    if (t->i < 0) {
                        throw std::runtime_error(std::string("Unknown variable ")
                                                 + ((struct parser_symbol*)node->r)->name);
                    }
                }
                t->reverse = true;
            }
            exe_size += sizeof(ParserExeDIV_PN);
        }
        else
        {
            int d1 = parser_ast_depth(node->l);
            int d2 = parser_ast_depth(node->r);
            if (d1 < d2) {
                parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                        local_variables);
                parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                        local_variables);
                if (p) {
                    new(p)      ParserExeDIV_B;
                    p += sizeof(ParserExeDIV_B);
                }
                exe_size += sizeof(ParserExeDIV_B);
            } else {
                parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                        local_variables);
                parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                        local_variables);
                if (p) {
                    new(p)      ParserExeDIV_F;
                    p += sizeof(ParserExeDIV_F);
                }
                exe_size += sizeof(ParserExeDIV_F);
            }
            --stack_size;
        }
        // no need to update max_grid_size
        break;
    }
    case PARSER_NEG:
    {
        parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size, local_variables);
        if (p) {
            new(p)      ParserExeNEG;
            p += sizeof(ParserExeNEG);
        }
        exe_size += sizeof(ParserExeNEG);
        break;
    }
    case PARSER_F1:
    {
        parser_compile_exe_size(((struct parser_f1*)node)->l,
                                p, exe_size, max_stack_size, stack_size, local_variables);
        if (p) {
            auto *t = new(p) ParserExeF1;
            p     += sizeof(ParserExeF1);
            t->ftype = ((struct parser_f1*)node)->ftype;
        }
        exe_size += sizeof(ParserExeF1);
        break;
    }
    case PARSER_F2:
    {
        int d1 = parser_ast_depth(((struct parser_f2*)node)->l);
        int d2 = parser_ast_depth(((struct parser_f2*)node)->r);
        if (d1 < d2) {
            parser_compile_exe_size(((struct parser_f2*)node)->r,
                                    p, exe_size, max_stack_size, stack_size, local_variables);
            parser_compile_exe_size(((struct parser_f2*)node)->l,
                                    p, exe_size, max_stack_size, stack_size, local_variables);
            if (p) {
                auto *t = new(p) ParserExeF2_B;
                p     += sizeof(ParserExeF2_B);
                t->ftype = ((struct parser_f2*)node)->ftype;
            }
            exe_size += sizeof(ParserExeF2_B);
        } else {
            parser_compile_exe_size(((struct parser_f2*)node)->l,
                                    p, exe_size, max_stack_size, stack_size, local_variables);
            parser_compile_exe_size(((struct parser_f2*)node)->r,
                                    p, exe_size, max_stack_size, stack_size, local_variables);
            if (p) {
                auto *t = new(p) ParserExeF2_F;
                p     += sizeof(ParserExeF2_F);
                t->ftype = ((struct parser_f2*)node)->ftype;
            }
            exe_size += sizeof(ParserExeF2_F);
        }
        --stack_size;
        break;
    }
    case PARSER_F3:
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(((struct parser_f3*)node)->ftype == PARSER_IF,
                                         "parser_compile: unknown f3 type");
        parser_compile_exe_size(((struct parser_f3*)node)->n1,
                                p, exe_size, max_stack_size, stack_size, local_variables);

        ParserExeIF* tif = nullptr;
        char* psave = nullptr;
        if (p) {
            tif = new(p) ParserExeIF;
            p  += sizeof(ParserExeIF);
            psave = p;
        }
        exe_size += sizeof(ParserExeIF);
        --stack_size;
        auto stack_size_save = stack_size;

        parser_compile_exe_size(((struct parser_f3*)node)->n2,
                                p, exe_size, max_stack_size, stack_size, local_variables);

        ParserExeJUMP* tjump = nullptr;
        if (p) {
            tjump = new(p) ParserExeJUMP;
            p +=    sizeof(ParserExeJUMP);
        }
        exe_size += sizeof(ParserExeJUMP);

        if (psave) {
            tif->offset = static_cast<int>(p-psave);
        }
        stack_size = stack_size_save;

        psave = p;
        parser_compile_exe_size(((struct parser_f3*)node)->n3,
                                p, exe_size, max_stack_size, stack_size, local_variables);
        if (tjump) {
            tjump->offset = static_cast<int>(p-psave);
        }

        break;
    }
    case PARSER_ASSIGN:
    {
        auto *asgn = (struct parser_assign*)node;
        local_variables.push_back(asgn->s->name);
        parser_compile_exe_size(asgn->v, p, exe_size, max_stack_size, stack_size, local_variables);
        break;
    }
    case PARSER_LIST:
    {
        parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size, local_variables);
        parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size, local_variables);
        break;
    }
    case PARSER_ADD_VP:
    {
        if (p) {
            auto *t = new(p) ParserExeADD_VP;
            p     += sizeof(ParserExeADD_VP);
            int lidx = parser_local_symbol_index((struct parser_symbol*)(node->r), local_variables);
            if (lidx >= 0) {
                t->i = AMREX_PARSER_LOCAL_IDX0 + lidx;
            } else {
                t->i = node->rip;
                if (t->i < 0) {
                  throw std::runtime_error(std::string("Unknown variable ")
                                           + ((struct parser_symbol*)node->r)->name);
                }
            }
            t->v = node->lvp.v;
        }
        exe_size += sizeof(ParserExeADD_VP);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case PARSER_SUB_VP:
    {
        if (p) {
            auto *t = new(p) ParserExeSUB_VP;
            p     += sizeof(ParserExeSUB_VP);
            int lidx = parser_local_symbol_index((struct parser_symbol*)(node->r), local_variables);
            if (lidx >= 0) {
                t->i = AMREX_PARSER_LOCAL_IDX0 + lidx;
            } else {
                t->i = node->rip;
                if (t->i < 0) {
                  throw std::runtime_error(std::string("Unknown variable ")
                                           + ((struct parser_symbol*)node->r)->name);
                }
            }
            t->v = node->lvp.v;
        }
        exe_size += sizeof(ParserExeSUB_VP);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case PARSER_MUL_VP:
    {
        if (p) {
            auto *t = new(p) ParserExeMUL_VP;
            p     += sizeof(ParserExeMUL_VP);
            int lidx = parser_local_symbol_index((struct parser_symbol*)(node->r), local_variables);
            if (lidx >= 0) {
                t->i = AMREX_PARSER_LOCAL_IDX0 + lidx;
            } else {
                t->i = node->rip;
                if (t->i < 0) {
                  throw std::runtime_error(std::string("Unknown variable ")
                                           + ((struct parser_symbol*)node->r)->name);
                }
            }
            t->v = node->lvp.v;
        }
        exe_size += sizeof(ParserExeMUL_VP);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case PARSER_DIV_VP:
    {
        if (p) {
            auto *t = new(p) ParserExeDIV_VP;
            p     += sizeof(ParserExeDIV_VP);
            int lidx = parser_local_symbol_index((struct parser_symbol*)(node->r), local_variables);
            if (lidx >= 0) {
                t->i = AMREX_PARSER_LOCAL_IDX0 + lidx;
            } else {
                t->i = node->rip;
                if (t->i < 0) {
                  throw std::runtime_error(std::string("Unknown variable ")
                                           + ((struct parser_symbol*)node->r)->name);
                }
            }
            t->v = node->lvp.v;
        }
        exe_size += sizeof(ParserExeDIV_VP);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case PARSER_ADD_PP:
    {
        if (p) {
            auto *t = new(p) ParserExeADD_PP;
            p     += sizeof(ParserExeADD_PP);
            int li1 = parser_local_symbol_index((struct parser_symbol*)(node->l), local_variables);
            int li2 = parser_local_symbol_index((struct parser_symbol*)(node->r), local_variables);
            if (li1 >= 0) {
                t->i1 = AMREX_PARSER_LOCAL_IDX0 + li1;
            } else {
                t->i1 = node->lvp.ip;
            }
            if (li2 >= 0) {
                t->i2 = AMREX_PARSER_LOCAL_IDX0 + li2;
            } else {
                t->i2 = node->rip;
            }
        }
        exe_size += sizeof(ParserExeADD_PP);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case PARSER_SUB_PP:
    {
        if (p) {
            auto *t = new(p) ParserExeSUB_PP;
            p     += sizeof(ParserExeSUB_PP);
            int li1 = parser_local_symbol_index((struct parser_symbol*)(node->l), local_variables);
            int li2 = parser_local_symbol_index((struct parser_symbol*)(node->r), local_variables);
            if (li1 >= 0) {
                t->i1 = AMREX_PARSER_LOCAL_IDX0 + li1;
            } else {
                t->i1 = node->lvp.ip;
            }
            if (li2 >= 0) {
                t->i2 = AMREX_PARSER_LOCAL_IDX0 + li2;
            } else {
                t->i2 = node->rip;
            }
        }
        exe_size += sizeof(ParserExeSUB_PP);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case PARSER_MUL_PP:
    {
        if (p) {
            auto *t = new(p) ParserExeMUL_PP;
            p     += sizeof(ParserExeMUL_PP);
            int li1 = parser_local_symbol_index((struct parser_symbol*)(node->l), local_variables);
            int li2 = parser_local_symbol_index((struct parser_symbol*)(node->r), local_variables);
            if (li1 >= 0) {
                t->i1 = AMREX_PARSER_LOCAL_IDX0 + li1;
            } else {
                t->i1 = node->lvp.ip;
            }
            if (li2 >= 0) {
                t->i2 = AMREX_PARSER_LOCAL_IDX0 + li2;
            } else {
                t->i2 = node->rip;
            }
        }
        exe_size += sizeof(ParserExeMUL_PP);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case PARSER_DIV_PP:
    {
        if (p) {
            auto *t = new(p) ParserExeDIV_PP;
            p     += sizeof(ParserExeDIV_PP);
            int li1 = parser_local_symbol_index((struct parser_symbol*)(node->l), local_variables);
            int li2 = parser_local_symbol_index((struct parser_symbol*)(node->r), local_variables);
            if (li1 >= 0) {
                t->i1 = AMREX_PARSER_LOCAL_IDX0 + li1;
            } else {
                t->i1 = node->lvp.ip;
            }
            if (li2 >= 0) {
                t->i2 = AMREX_PARSER_LOCAL_IDX0 + li2;
            } else {
                t->i2 = node->rip;
            }
        }
        exe_size += sizeof(ParserExeDIV_PP);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case PARSER_NEG_P:
    {
        if (p) {
            auto *t = new(p) ParserExeNEG_P;
            p     += sizeof(ParserExeNEG_P);
            int lidx = parser_local_symbol_index((struct parser_symbol*)(node->l), local_variables);
            if (lidx >= 0) {
                t->i = AMREX_PARSER_LOCAL_IDX0 + lidx;
            } else {
                t->i = node->lvp.ip;
                if (t->i < 0) {
                  throw std::runtime_error(std::string("Unknown variable ")
                                           + ((struct parser_symbol*)node->l)->name);
                }
            }
        }
        exe_size += sizeof(ParserExeNEG_P);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    default:
        amrex::Abort("parser_compile: unknown node type " + std::to_string(node->type));
    }
}

}
