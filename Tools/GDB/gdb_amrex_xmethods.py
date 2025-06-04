# pylint: disable=too-few-public-methods, too-many-arguments, too-many-boolean-expressions

import os.path
import sys

sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "util"))
import class_methods


def is_f_order(obj):
    return obj.ORDER.name == "amrex::Order::F"


def is_c_order(obj):
    return obj.ORDER.name == "amrex::Order::C"


@class_methods.Class("amrex::Array1D", template_types=["T", "XLO", "XHI"])
class AmrexArray1D:  # pylint: disable=no-member
    @class_methods.member_function("unsigned int", "len", [])
    def len(self, _obj):
        return self.XHI - self.XLO + 1

    @class_methods.member_function("unsigned int", "size", [])
    def size(self, obj):
        return self.len(obj)

    @class_methods.member_function("T&", "operator()", ["int"])
    def subscript(self, obj, i):
        if not self.XLO <= i <= self.XHI:
            msg = 'Array1D index "{}" should be between {} and {}.'.format(
                i, self.XLO, self.XHI
            )
            raise IndexError(msg)
        return obj["arr"][i - self.XLO]


@class_methods.Class(
    "amrex::Array2D", template_types=["T", "XLO", "XHI", "YLO", "YHI", "ORDER"]
)
class AmrexArray2D:  # pylint: disable=no-member
    @class_methods.member_function("unsigned int", "xlen", [])
    def xlen(self, _obj):
        return self.XHI - self.XLO + 1

    @class_methods.member_function("unsigned int", "ylen", [])
    def ylen(self, _obj):
        return self.YHI - self.YLO + 1

    @class_methods.member_function("unsigned int", "size", [])
    def size(self, obj):
        return self.xlen(obj) * self.ylen(obj)

    @class_methods.member_function("T&", "operator()", ["int", "int"])
    def subscript(self, obj, i, j):
        if not self.XLO <= i <= self.XHI:
            msg = 'Array2D index i="{}" should be between {} and {}.'.format(
                i, self.XLO, self.XHI
            )
            raise IndexError(msg)
        if not self.YLO <= j <= self.YHI:
            msg = 'Array2D index j="{}" should be between {} and {}.'.format(
                j, self.YLO, self.YHI
            )
            raise IndexError(msg)

        xlen = self.xlen(obj)
        ylen = self.ylen(obj)
        if is_f_order(self):
            return obj["arr"][i + j * xlen - (self.YLO * xlen + self.XLO)]
        if is_c_order(self):
            return obj["arr"][j + i * ylen - (self.XLO * ylen + self.YLO)]
        assert False


@class_methods.Class(
    "amrex::Array3D",
    template_types=["T", "XLO", "XHI", "YLO", "YHI", "ZLO", "ZHI", "ORDER"],
)
class AmrexArray3D:  # pylint: disable=no-member
    @class_methods.member_function("unsigned int", "xlen", [])
    def xlen(self, _obj):
        return self.XHI - self.XLO + 1

    @class_methods.member_function("unsigned int", "ylen", [])
    def ylen(self, _obj):
        return self.YHI - self.YLO + 1

    @class_methods.member_function("unsigned int", "zlen", [])
    def zlen(self, _obj):
        return self.ZHI - self.ZLO + 1

    @class_methods.member_function("unsigned int", "size", [])
    def size(self, obj):
        return self.xlen(obj) * self.ylen(obj) * self.zlen(obj)

    @class_methods.member_function("T&", "operator()", ["int", "int", "int"])
    def element(self, obj, i, j, k):
        if not self.XLO <= i <= self.XHI:
            msg = 'Array3D index i="{}" should be between {} and {}.'.format(
                i, self.XLO, self.XHI
            )
            raise IndexError(msg)
        if not self.YLO <= j <= self.YHI:
            msg = 'Array3D index j="{}" should be between {} and {}.'.format(
                j, self.YLO, self.YHI
            )
            raise IndexError(msg)
        if not self.ZLO <= k <= self.ZHI:
            msg = 'Array3D index k="{}" should be between {} and {}.'.format(
                k, self.ZLO, self.ZHI
            )
            raise IndexError(msg)

        xlen = self.xlen(obj)
        ylen = self.ylen(obj)
        zlen = self.zlen(obj)
        if is_f_order(self):
            return obj["arr"][
                i
                + j * xlen
                + k * (xlen * ylen)
                - (self.ZLO * (xlen * ylen) + self.YLO * xlen + self.XLO)
            ]
        if is_c_order(self):
            return obj["arr"][
                k
                + j * zlen
                + i * (zlen * ylen)
                - (self.XLO * (zlen * ylen) + self.YLO * zlen + self.ZLO)
            ]
        assert False


@class_methods.Class("amrex::Array4", template_types=["T"])
class AmrexArray4:
    @class_methods.member_function("std::size_t", "size", [])
    def size(self, obj):
        return obj["nstride"] * obj["ncomp"]

    def subscript_helper(self, obj, i, j, k, n):
        begin = obj["begin"]
        end = obj["end"]
        # fmt: off
        if (
            i < begin["x"] or i >= end["x"] or
            j < begin["y"] or j >= end["y"] or
            k < begin["z"] or k >= end["z"] or
            n < 0 or n >= obj["ncomp"]
        ):
            msg = "Array4 index ({},{},{},{}) is out of bound ({}:{},{}:{},{}:{},0:{})".format(
                i, j, k, n,
                begin["x"], end["x"] - 1,
                begin["y"], end["y"] - 1,
                begin["z"], end["z"] - 1,
                obj["ncomp"] - 1,
            )
            raise IndexError(msg)
        # fmt: on
        return obj["p"][
            (i - begin["x"])
            + (j - begin["y"]) * obj["jstride"]
            + (k - begin["z"]) * obj["kstride"]
            + n * obj["nstride"]
        ]

    # GDB's overload resolution ignores any missing arguments at the end, so
    # this handles both the 3-int and 4-int overloads
    @class_methods.member_function("T&", "operator()", ["int"] * 4)
    def subscript(self, obj, i, j, k, n=0):
        return self.subscript_helper(obj, i, j, k, n)

    @class_methods.member_function("T&", "operator()", ["amrex::IntVect&", "int"])
    def subscript_IntVect(self, obj, iv, n=0):
        spacedim = int(iv.type.template_argument(0))
        indices = [iv["vect"][0], 0, 0]
        if spacedim >= 2:
            indices[1] = iv["vect"][1]
        if spacedim >= 3:
            indices[2] = iv["vect"][2]
        return self.subscript_helper(obj, *indices, n)

    @class_methods.member_function("T&", "operator()", ["amrex::Dim3&", "int"])
    def subscript_Dim3(self, obj, cell, n=0):
        return self.subscript_helper(obj, cell["x"], cell["y"], cell["z"], n)


@class_methods.Class("amrex::Table1D", template_types=["T"])
class AmrexTable1D:
    @class_methods.member_function("T&", "operator()", ["int"])
    def subscript(self, obj, i):
        begin = obj["begin"]
        end = obj["end"]
        if i < begin or i >= end:
            msg = f"({i}) is out of bound ({begin}:{end-1})"
            raise IndexError(msg)
        return obj["p"][i - begin]


@class_methods.Class("amrex::Table2D", template_types=["T", "ORDER"])
class AmrexTable2D:
    @class_methods.member_function("T&", "operator()", ["int"] * 2)
    def subscript(self, obj, i, j):
        begin = obj["begin"]
        end = obj["end"]
        # fmt: off
        if (
            i < begin[0] or i >= end[0] or
            j < begin[1] or j >= end[1]
        ):
            msg = "({},{}) is out of bound ({}:{},{}:{})".format(
                i, j,
                begin[0], end[0] - 1,
                begin[1], end[1] - 1,
            )
            raise IndexError(msg)
        # fmt: on

        stride1 = obj["stride1"]
        if is_f_order(self):
            return obj["p"][(i - begin[0]) + (j - begin[1]) * stride1]
        if is_c_order(self):
            return obj["p"][(i - begin[0]) * stride1 + (j - begin[1])]
        assert False


@class_methods.Class("amrex::Table3D", template_types=["T", "ORDER"])
class AmrexTable3D:
    @class_methods.member_function("T&", "operator()", ["int"] * 3)
    def subscript(self, obj, i, j, k):
        begin = obj["begin"]["arr"]
        end = obj["end"]["arr"]
        # fmt: off
        if (
            i < begin[0] or i >= end[0] or
            j < begin[1] or j >= end[1] or
            k < begin[2] or k >= end[2]
        ):
            msg = "({},{},{}) is out of bound ({}:{},{}:{},{}:{})".format(
                i, j, k,
                begin[0], end[0] - 1,
                begin[1], end[1] - 1,
                begin[2], end[2] - 1,
            )
            raise ValueError(msg)
        # fmt: on

        stride1 = obj["stride1"]
        stride2 = obj["stride2"]
        if is_f_order(self):
            return obj["p"][
                (i - begin[0]) + (j - begin[1]) * stride1 + (k - begin[2]) * stride2
            ]
        if is_c_order(self):
            return obj["p"][
                (i - begin[0]) * stride2 + (j - begin[1]) * stride1 + (k - begin[2])
            ]
        assert False


@class_methods.Class("amrex::Table4D", template_types=["T", "ORDER"])
class AmrexTable4D:
    @class_methods.member_function("T&", "operator()", ["int"] * 4)
    def subscript(self, obj, i, j, k, n):
        begin = obj["begin"]["arr"]
        end = obj["end"]["arr"]
        # fmt: off
        if (
            i < begin[0] or i >= end[0] or
            j < begin[1] or j >= end[1] or
            k < begin[2] or k >= end[2] or
            n < begin[3] or n >= end[3]
        ):
            msg = "({},{},{},{}) is out of bound ({}:{},{}:{},{}:{},{}:{})".format(
                i, j, k, n,
                begin[0], end[0] - 1,
                begin[1], end[1] - 1,
                begin[2], end[2] - 1,
                begin[3], end[3] - 1,
            )
            raise ValueError(msg)
        # fmt: on

        stride1 = obj["stride1"]
        stride2 = obj["stride2"]
        stride3 = obj["stride3"]
        if is_f_order(self):
            return obj["p"][
                (i - begin[0])
                + (j - begin[1]) * stride1
                + (k - begin[2]) * stride2
                + (n - begin[3]) * stride3
            ]
        if is_c_order(self):
            return obj["p"][
                (i - begin[0]) * stride3
                + (j - begin[1]) * stride2
                + (k - begin[2]) * stride1
                + (n - begin[3])
            ]
        assert False
