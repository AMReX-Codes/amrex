! MODULE mod_interpret
! Parser to read and evaluate math expressions

module mod_interpret
use iso_c_binding
use amrex_fort_module, only : amrex_real
implicit none
INTEGER, parameter :: plus       =1, &
                      minus      =2, &
                      multiply   =3, &
                      divide     =4, &
                      power      =5, &
                      greaterthan=6, &
                      lessthan   =7, &
                      exponential=8, &
                      logarithm  =9, &
                      sine       =10, &
                      cosine     =11, &
                      tangent    =12, &
                      square_root=13, &
                      arccos     =14, &
                      arcsin     =15, &
                      arctan     =16, &
                      sinhyp     =17, &
                      coshyp     =18, &
                      tanhyp     =19, &
                      logten     =20

CHARACTER(5), DIMENSION(0:20) :: coper
INTEGER, parameter :: w_parenthesis = 1, &
                      w_number      = 2, &
                      w_operator    = 3, &
                      w_intrinsic   = 4, &
                      w_variable    = 5
TYPE :: operation_type
  INTEGER :: a,b,op
END TYPE operation_type

TYPE :: res_type
  INTEGER :: nb_op
  REAL(amrex_real) :: value
  TYPE(res_type), DIMENSION(:), POINTER :: res
  LOGICAL, DIMENSION(:), POINTER :: l_res, l_res_do
  TYPE(operation_type), DIMENSION(:), POINTER :: operation
  LOGICAL :: a_res,a_l_res,a_l_res_do,a_operation
END TYPE res_type

TYPE :: res_type_array
  INTEGER :: nb_op
  REAL(amrex_real), DIMENSION(:), POINTER :: value
  TYPE(res_type_array), DIMENSION(:), POINTER :: res
  LOGICAL, DIMENSION(:), POINTER :: l_res, l_res_do
  TYPE(operation_type), DIMENSION(:), POINTER :: operation
  LOGICAL :: a_res,a_l_res,a_l_res_do,a_operation
END TYPE res_type_array

TYPE(res_type_array) :: res_array

INTEGER :: nb_res=0
TYPE(res_type), DIMENSION(10) :: table_of_res

contains

RECURSIVE subroutine del_res(resin)
implicit none
TYPE(res_type), INTENT(INOUT) :: resin
INTEGER :: j
  IF(resin%a_res) then
    do j = LBOUND(resin%res,1),UBOUND(resin%res,1)
      call del_res(resin%res(j))
    end do
   DEALLOCATE(resin%res)
   resin%a_res=.FALSE.
  END if
  IF(resin%a_l_res) then
    DEALLOCATE(resin%l_res)
    resin%a_l_res=.FALSE.
  END if
  IF(resin%a_l_res_do) then
    DEALLOCATE(resin%l_res_do)
    resin%a_l_res_do=.FALSE.
  END if
  IF(resin%a_operation) then
    DEALLOCATE(resin%operation)
    resin%a_operation=.FALSE.
  END if
return
end subroutine del_res

RECURSIVE subroutine del_res_array(resin)
implicit none
TYPE(res_type_array), INTENT(INOUT) :: resin
INTEGER :: j
  IF(resin%a_res) then
    do j = LBOUND(resin%res,1),UBOUND(resin%res,1)
      call del_res_array(resin%res(j))
    end do
   DEALLOCATE(resin%res)
   resin%a_res=.FALSE.
  END if
  IF(resin%a_l_res) then
    DEALLOCATE(resin%l_res)
    resin%a_l_res=.FALSE.
  END if
  IF(resin%a_l_res_do) then
    DEALLOCATE(resin%l_res_do)
    resin%a_l_res_do=.FALSE.
  END if
  IF(resin%a_operation) then
    DEALLOCATE(resin%operation)
    resin%a_operation=.FALSE.
  END if
  if (associated(resin%value)) then
      nullify(resin%value)
  end if
return
end subroutine del_res_array

subroutine init_interpret()
implicit none

coper(0          )='###'
coper(plus       )='+'
coper(minus      )='-'
coper(multiply   )='*'
coper(divide     )='/'
coper(power      )='**'
coper(square_root)='sqrt'
coper(exponential)='exp'
coper(logarithm  )='log'
coper(sine       )='sin'
coper(cosine     )='cos'
coper(tangent    )='tan'
coper(arccos     )='acos'
coper(arcsin     )='asin'
coper(arctan     )='atan'
coper(sinhyp     )='sinh'
coper(coshyp     )='cosh'
coper(tanhyp     )='tanh'
coper(logten     )='log10'
coper(greaterthan)='>'
coper(lessthan   )='<'

return
end subroutine init_interpret

recursive function calc_res(res,list_var,root) RESULT(calc_res_r)
implicit none
REAL(amrex_real) :: calc_res_r
TYPE(res_type), INTENT(IN OUT) :: res
REAL(amrex_real), DIMENSION(:), OPTIONAL :: list_var
LOGICAL, INTENT(IN), OPTIONAL :: root

INTEGER :: j,ia,ib,op
REAL(amrex_real) :: a,b

do j = 1, res%nb_op
  ia = res%operation(j)%a
  ib = res%operation(j)%b
  op = res%operation(j)%op
  IF(op<0) then
    IF(op>-100) res%res(ia)%value=list_var(-op)
  else
    IF(res%l_res(ia).and.res%l_res_do(ia)) then
      IF(PRESENT(list_var)) then
        a = calc_res(res=res%res(ia),list_var=list_var,root=.false.)
      else
        a = calc_res(res=res%res(ia),root=.false.)
      END if
      res%l_res_do(ia)=.false.
    else
      a = res%res(ia)%value
    END if
    IF(res%l_res(ib).and.res%l_res_do(ib)) then
      IF(PRESENT(list_var)) then
        b = calc_res(res=res%res(ib),list_var=list_var,root=.false.)
      else
        b = calc_res(res=res%res(ib),root=.false.)
      END if
      ! This might prevent OMP threading
      res%l_res_do(ib)=.false.
    else
      b = res%res(ib)%value
    END if
    ! This might prevent OMP threading
    res%res(ia)%value=eval(a,b,op)
  END if
end do
IF(res%nb_op==0) then
  ia=1
  calc_res_r = res%res(ia)%value
else
  calc_res_r = res%res(ia)%value
END if

IF(.not.PRESENT(root))  call calc_res_init(res)

end function calc_res

recursive subroutine calc_res_init(res)
implicit none
TYPE(res_type), INTENT(IN OUT) :: res

INTEGER :: j,ia,ib,op

do j = 1, res%nb_op
  ia = res%operation(j)%a
  ib = res%operation(j)%b
  op = res%operation(j)%op
  IF(op>=0) then
    IF(res%l_res(ia).and..not.res%l_res_do(ia)) then
      call calc_res_init(res%res(ia))
      res%l_res_do(ia)=.true.
    END if
    IF(res%l_res(ib).and..not.res%l_res_do(ib)) then
      call calc_res_init(res%res(ib))
      res%l_res_do(ib)=.true.
    END if
  END if
end do

end subroutine calc_res_init

recursive function calc_res_array(res,list_values,root) RESULT(calc_res_array_r)
implicit none
REAL(amrex_real), dimension(:), pointer :: calc_res_array_r
TYPE(res_type_array), INTENT(IN OUT) :: res
REAL(amrex_real), DIMENSION(:,:), OPTIONAL :: list_values
LOGICAL, INTENT(IN), OPTIONAL :: root

INTEGER :: j,ia,ib,op,n
REAL(amrex_real),dimension(:),pointer :: a,b

do j = 1, res%nb_op
  ia = res%operation(j)%a
  ib = res%operation(j)%b
  op = res%operation(j)%op
  IF(op<0) then
    IF(op>-100) then
       n = size(list_values,2)
       allocate(res%res(ia)%value(n))
       res%res(ia)%value=list_values(-op,:)
    end if  
  else
    IF(res%l_res(ia).and.res%l_res_do(ia)) then
      IF(PRESENT(list_values)) then
        a => calc_res_array(res=res%res(ia),list_values=list_values,root=.false.)
      else
        a => calc_res_array(res=res%res(ia),root=.false.)
      END if
      res%l_res_do(ia)=.false.
    else
      a => res%res(ia)%value
    END if
    IF(res%l_res(ib).and.res%l_res_do(ib)) then
      IF(PRESENT(list_values)) then
        b => calc_res_array(res=res%res(ib),list_values=list_values,root=.false.)
      else
        b => calc_res_array(res=res%res(ib),root=.false.)
      END if
      res%l_res_do(ib)=.false.
    else
      b => res%res(ib)%value
    END if
    res%res(ia)%value=>eval_array(a,b,op)
  END if
end do

IF(res%nb_op==0) then
  ia=1
  calc_res_array_r => res%res(ia)%value
else
  calc_res_array_r => res%res(ia)%value
END if

IF(.not.PRESENT(root))  call calc_res_array_init(res)

end function calc_res_array

recursive subroutine calc_res_array_init(res)
implicit none
TYPE(res_type_array), INTENT(IN OUT) :: res

INTEGER :: j,ia,ib,op

do j = 1, res%nb_op
  ia = res%operation(j)%a
  ib = res%operation(j)%b
  op = res%operation(j)%op
  IF(op>=0) then
    IF(res%l_res(ia).and..not.res%l_res_do(ia)) then
      call calc_res_array_init(res%res(ia))
      res%l_res_do(ia)=.true.
    END if
    IF(res%l_res(ib).and..not.res%l_res_do(ib)) then
      call calc_res_array_init(res%res(ib))
      res%l_res_do(ib)=.true.
    END if
  END if
end do

end subroutine calc_res_array_init

recursive subroutine eval_res(res,exprin,int_op,list_var,l_root)
implicit none
TYPE(res_type) :: res
CHARACTER(LEN=*), INTENT(IN) :: exprin
INTEGER, OPTIONAL, INTENT(IN) :: int_op
CHARACTER(*), DIMENSION(:),INTENT(IN), OPTIONAL :: list_var
LOGICAL, OPTIONAL :: l_root

CHARACTER(LEN=LEN(exprin)) :: expr_old
CHARACTER(LEN=LEN(exprin)+2) :: expr
INTEGER, parameter :: nop_tot=25
INTEGER :: ln, i, j,ja
INTEGER :: jtot,jres,jop
INTEGER, DIMENSION(0:nop_tot) :: next_l,next_r,op,what
LOGICAL, ALLOCATABLE, DIMENSION(:) :: l_res

IF(.not.PRESENT(l_root)) call del_res(res)

res%a_res=.false.
res%a_l_res=.false.
res%a_l_res_do=.false.
res%a_operation=.false.

op=0
what=0

expr_old = TRIM(ADJUSTL(C_replace(exprin,' ','')))
IF(expr(1:1)=='-') then
  expr='0'//expr_old
else
  expr='0+'//expr_old
END if
ln = LEN(TRIM(ADJUSTL(expr)))
next_r(0) = 0
i = 1
j = 1
jres = 0
do WHILE(i<ln+1)
  IF(PRESENT(list_var)) then
    call what_is_next(expr(i:ln),what(j),next_l(j),next_r(j),op(j),list_var)
  else
    call what_is_next(expr(i:ln),what(j),next_l(j),next_r(j),op(j))
  END if
  next_r(j) = next_r(j)+next_r(j-1)
  next_l(j) = next_r(j-1)+1+next_l(j)
  i = next_r(j)+1
  IF(what(j)/=w_operator) jres=jres+1
  j = j+1
end do
jtot = j-1
IF(PRESENT(int_op)) then
  ALLOCATE(res%res(0:jres),res%operation(jres),res%l_res(0:jres),res%l_res_do(0:jres))
  res%a_res=.true.
  res%a_l_res=.true.
  res%a_l_res_do=.true.
  res%a_operation=.true.
  res%res(0:jres)%a_res=.false.
  res%res(0:jres)%a_l_res=.false.
  res%res(0:jres)%a_l_res_do=.false.
  res%res(0:jres)%a_operation=.false.
  res%l_res(0)=.false.
  res%l_res_do=.true.
  res%res(0)%value=0.
  res%nb_op=jres
  res%operation(jres)%op=int_op
  res%operation(jres)%a=1
  res%operation(jres)%b=0
  res%l_res(jres)=.false.
  nullify(res%res(0)%res)
  nullify(res%res(0)%l_res)
  nullify(res%res(0)%l_res_do)
  nullify(res%res(0)%operation)
else
  ALLOCATE(res%res(jres),res%operation(jres-1),res%l_res(jres),res%l_res_do(jres))
  res%a_res=.true.
  res%a_l_res=.true.
  res%a_l_res_do=.true.
  res%a_operation=.true.
  res%res(1:jres)%a_res=.false.
  res%res(1:jres)%a_l_res=.false.
  res%res(1:jres)%a_l_res_do=.false.
  res%res(1:jres)%a_operation=.false.
  res%nb_op=jres-1
  res%l_res_do=.true.
END if
ALLOCATE(l_res(jres))
l_res(:) = .true.

res%res(:)%value=0.

jres = 0
do j=1, jtot
  select case (what(j))
    case (w_number)
      jres = jres+1
      ALLOCATE(res%res(jres)%res(1),res%res(jres)%l_res(1),res%res(jres)%l_res_do(1),res%res(jres)%operation(1))
      res%res(jres)%a_res=.true.
      res%res(jres)%a_l_res=.true.
      res%res(jres)%a_l_res_do=.true.
      res%res(jres)%a_operation=.true.
      res%res(jres)%res(1)%a_res=.false.
      res%res(jres)%res(1)%a_l_res=.false.
      res%res(jres)%res(1)%a_l_res_do=.false.
      res%res(jres)%res(1)%a_operation=.false.
      res%res(jres)%nb_op=1
      res%res(jres)%l_res(1) = .false.
      res%res(jres)%l_res_do(1) = .true.
      res%res(jres)%operation(1)%op = -100-op(j)
      res%res(jres)%operation(1)%a = 1
      res%res(jres)%operation(1)%b = 0
      res%l_res(jres)=.true.
      res%res(jres)%res(1)%value = char2real(expr(next_l(j):next_r(j)))
      nullify(res%res(jres)%res(1)%res)
      nullify(res%res(jres)%res(1)%l_res)
      nullify(res%res(jres)%res(1)%l_res_do)
      nullify(res%res(jres)%res(1)%operation)
    case (w_parenthesis)
      jres = jres+1
      IF(PRESENT(list_var)) then
        call eval_res(res%res(jres),expr(next_l(j)+1:next_r(j)-1),list_var=list_var,l_root=.false.)
      else
        call eval_res(res%res(jres),expr(next_l(j)+1:next_r(j)-1),l_root=.false.)
      endif
      res%l_res(jres)=.true.
    case (w_intrinsic)
      jres = jres+1
      IF(PRESENT(list_var)) then
        call eval_res(res%res(jres),expr(next_l(j)+1:next_r(j)-1),op(j),list_var,l_root=.false.)
      else
        call eval_res(res%res(jres),expr(next_l(j)+1:next_r(j)-1),op(j),l_root=.false.)
      END if
      res%l_res(jres)=.true.
    case (w_variable)
      jres = jres+1
      ALLOCATE(res%res(jres)%res(1),res%res(jres)%l_res(1),res%res(jres)%l_res_do(1),res%res(jres)%operation(1))
      res%res(jres)%a_res=.true.
      res%res(jres)%a_l_res=.true.
      res%res(jres)%a_l_res_do=.true.
      res%res(jres)%a_operation=.true.
      res%res(jres)%res(1)%a_res=.false.
      res%res(jres)%res(1)%a_l_res=.false.
      res%res(jres)%res(1)%a_l_res_do=.false.
      res%res(jres)%res(1)%a_operation=.false.
      res%res(jres)%nb_op=1
      res%res(jres)%l_res(1) = .false.
      res%res(jres)%l_res_do(1) = .true.
      res%res(jres)%operation(1)%op = -op(j)
      res%res(jres)%operation(1)%a = 1
      res%res(jres)%operation(1)%b = 0
      res%l_res(jres)=.true.
      nullify(res%res(jres)%res(1)%res)
      nullify(res%res(jres)%res(1)%l_res)
      nullify(res%res(jres)%res(1)%l_res_do)
      nullify(res%res(jres)%res(1)%operation)
    case (w_operator)
    case default
  end select
end do
jop = 1
do j=jtot,1,-1
  IF(what(j)==w_operator) then
    IF(op(j)==power) then
      res%operation(jop)%op = power
      ja=j/2
      do WHILE(.not.l_res(ja))
        ja=ja-1
      end do
      res%operation(jop)%a = ja
      res%operation(jop)%b = j/2+1
      l_res(j/2+1)=.false.
      jop = jop+1
    END IF
  END if
end do
do j=jtot,1,-1
  IF(what(j)==w_operator) then
    IF(op(j)==multiply.or.op(j)==divide) then
      res%operation(jop)%op = op(j)
      ja=j/2
      do WHILE(.not.l_res(ja))
        ja=ja-1
      end do
      res%operation(jop)%a = ja
      res%operation(jop)%b = j/2+1
      l_res(j/2+1)=.false.
      jop = jop+1
    END if
  END if
end do
do j=jtot,1,-1
  IF(what(j)==w_operator) then
    IF(op(j)==plus.or.op(j)==minus) then
      res%operation(jop)%op = op(j)
      ja=j/2
      do WHILE(.not.l_res(ja))
        ja=ja-1
      end do
      res%operation(jop)%a = ja
      res%operation(jop)%b = j/2+1
      l_res(j/2+1)=.false.
      jop = jop+1
    END if
  END if
end do
do j=jtot,1,-1
  IF(what(j)==w_operator) then
    IF(op(j)==greaterthan.or.op(j)==lessthan) then
      res%operation(jop)%op = op(j)
      ja=j/2
      do WHILE(.not.l_res(ja))
        ja=ja-1
      end do
      res%operation(jop)%a = ja
      res%operation(jop)%b = j/2+1
      l_res(j/2+1)=.false.
      jop = jop+1
    END if
  END if
end do

DEALLOCATE(l_res)

return
end subroutine eval_res

recursive subroutine eval_res_array(res,exprin,int_op,list_var,l_root)
implicit none
TYPE(res_type_array) :: res
CHARACTER(LEN=*), INTENT(IN) :: exprin
INTEGER, OPTIONAL, INTENT(IN) :: int_op
CHARACTER(*), DIMENSION(:),INTENT(IN), OPTIONAL :: list_var
LOGICAL, OPTIONAL :: l_root

CHARACTER(LEN=LEN(exprin)) :: expr_old
CHARACTER(LEN=LEN(exprin)+2) :: expr
INTEGER, parameter :: nop_tot=25
INTEGER :: ln, i, j,ja
INTEGER :: jtot,jres,jop
INTEGER, DIMENSION(0:nop_tot) :: next_l,next_r,op,what
LOGICAL, ALLOCATABLE, DIMENSION(:) :: l_res

IF(.not.PRESENT(l_root)) call del_res_array(res)

res%a_res=.false.
res%a_l_res=.false.
res%a_l_res_do=.false.
res%a_operation=.false.

op=0
what=0

expr_old = TRIM(ADJUSTL(C_replace(exprin,' ','')))
IF(expr(1:1)=='-') then
  expr='0'//expr_old
else
  expr='0+'//expr_old
END if
ln = LEN(TRIM(ADJUSTL(expr)))
next_r(0) = 0
i = 1
j = 1
jres = 0
do WHILE(i<ln+1)
  IF(PRESENT(list_var)) then
    call what_is_next(expr(i:ln),what(j),next_l(j),next_r(j),op(j),list_var)
  else
    call what_is_next(expr(i:ln),what(j),next_l(j),next_r(j),op(j))
  END if
  next_r(j) = next_r(j)+next_r(j-1)
  next_l(j) = next_r(j-1)+1+next_l(j)
  i = next_r(j)+1
  IF(what(j)/=w_operator) jres=jres+1
  j = j+1
end do
jtot = j-1
IF(PRESENT(int_op)) then
  ALLOCATE(res%res(0:jres),res%operation(jres),res%l_res(0:jres),res%l_res_do(0:jres))
  res%a_res=.true.
  res%a_l_res=.true.
  res%a_l_res_do=.true.
  res%a_operation=.true.
  res%res(0:jres)%a_res=.false.
  res%res(0:jres)%a_l_res=.false.
  res%res(0:jres)%a_l_res_do=.false.
  res%res(0:jres)%a_operation=.false.
  res%l_res(0)=.false.
  res%l_res_do=.true.
  nullify(res%res(0)%value)
  res%nb_op=jres
  res%operation(jres)%op=int_op
  res%operation(jres)%a=1
  res%operation(jres)%b=0
  res%l_res(jres)=.false.
  nullify(res%res(0)%res)
  nullify(res%res(0)%l_res)
  nullify(res%res(0)%l_res_do)
  nullify(res%res(0)%operation)
else
  ALLOCATE(res%res(jres),res%operation(jres-1),res%l_res(jres),res%l_res_do(jres))
  res%a_res=.true.
  res%a_l_res=.true.
  res%a_l_res_do=.true.
  res%a_operation=.true.
  res%res(1:jres)%a_res=.false.
  res%res(1:jres)%a_l_res=.false.
  res%res(1:jres)%a_l_res_do=.false.
  res%res(1:jres)%a_operation=.false.
  res%nb_op=jres-1
  res%l_res_do=.true.
END if
ALLOCATE(l_res(jres))
l_res(:) = .true.

jres = 0
do j=1, jtot
  select case (what(j))
    case (w_number)
      jres = jres+1
      ALLOCATE(res%res(jres)%res(1),res%res(jres)%l_res(1),res%res(jres)%l_res_do(1),res%res(jres)%operation(1))
      res%res(jres)%a_res=.true.
      res%res(jres)%a_l_res=.true.
      res%res(jres)%a_l_res_do=.true.
      res%res(jres)%a_operation=.true.
      res%res(jres)%res(1)%a_res=.false.
      res%res(jres)%res(1)%a_l_res=.false.
      res%res(jres)%res(1)%a_l_res_do=.false.
      res%res(jres)%res(1)%a_operation=.false.
      res%res(jres)%nb_op=1
      res%res(jres)%l_res(1) = .false.
      res%res(jres)%l_res_do(1) = .true.
      res%res(jres)%operation(1)%op = -100-op(j)
      res%res(jres)%operation(1)%a = 1
      res%res(jres)%operation(1)%b = 0
      res%l_res(jres)=.true.
      allocate(res%res(jres)%res(1)%value(1))
      res%res(jres)%res(1)%value = char2real(expr(next_l(j):next_r(j)))
      nullify(res%res(jres)%res(1)%res)
      nullify(res%res(jres)%res(1)%l_res)
      nullify(res%res(jres)%res(1)%l_res_do)
      nullify(res%res(jres)%res(1)%operation)
    case (w_parenthesis)
      jres = jres+1
      IF(PRESENT(list_var)) then
        call eval_res_array(res%res(jres),expr(next_l(j)+1:next_r(j)-1),list_var=list_var,l_root=.false.)
      else
        call eval_res_array(res%res(jres),expr(next_l(j)+1:next_r(j)-1),l_root=.false.)
      endif
      res%l_res(jres)=.true.
    case (w_intrinsic)
      jres = jres+1
      IF(PRESENT(list_var)) then
        call eval_res_array(res%res(jres),expr(next_l(j)+1:next_r(j)-1),op(j),list_var,l_root=.false.)
      else
        call eval_res_array(res%res(jres),expr(next_l(j)+1:next_r(j)-1),op(j),l_root=.false.)
      END if
      res%l_res(jres)=.true.
    case (w_variable)
      jres = jres+1
      ALLOCATE(res%res(jres)%res(1),res%res(jres)%l_res(1),res%res(jres)%l_res_do(1),res%res(jres)%operation(1))
      res%res(jres)%a_res=.true.
      res%res(jres)%a_l_res=.true.
      res%res(jres)%a_l_res_do=.true.
      res%res(jres)%a_operation=.true.
      res%res(jres)%res(1)%a_res=.false.
      res%res(jres)%res(1)%a_l_res=.false.
      res%res(jres)%res(1)%a_l_res_do=.false.
      res%res(jres)%res(1)%a_operation=.false.
      res%res(jres)%nb_op=1
      res%res(jres)%l_res(1) = .false.
      res%res(jres)%l_res_do(1) = .true.
      res%res(jres)%operation(1)%op = -op(j)
      res%res(jres)%operation(1)%a = 1
      res%res(jres)%operation(1)%b = 0
      res%l_res(jres)=.true.
      nullify(res%res(jres)%res(1)%res)
      nullify(res%res(jres)%res(1)%l_res)
      nullify(res%res(jres)%res(1)%l_res_do)
      nullify(res%res(jres)%res(1)%operation)
    case (w_operator)
    case default
  end select
end do
jop = 1
do j=jtot,1,-1
  IF(what(j)==w_operator) then
    IF(op(j)==power) then
      res%operation(jop)%op = power
      ja=j/2
      do WHILE(.not.l_res(ja))
        ja=ja-1
      end do
      res%operation(jop)%a = ja
      res%operation(jop)%b = j/2+1
      l_res(j/2+1)=.false.
      jop = jop+1
    END IF
  END if
end do
do j=jtot,1,-1
  IF(what(j)==w_operator) then
    IF(op(j)==multiply.or.op(j)==divide) then
      res%operation(jop)%op = op(j)
      ja=j/2
      do WHILE(.not.l_res(ja))
        ja=ja-1
      end do
      res%operation(jop)%a = ja
      res%operation(jop)%b = j/2+1
      l_res(j/2+1)=.false.
      jop = jop+1
    END if
  END if
end do
do j=jtot,1,-1
  IF(what(j)==w_operator) then
    IF(op(j)==plus.or.op(j)==minus) then
      res%operation(jop)%op = op(j)
      ja=j/2
      do WHILE(.not.l_res(ja))
        ja=ja-1
      end do
      res%operation(jop)%a = ja
      res%operation(jop)%b = j/2+1
      l_res(j/2+1)=.false.
      jop = jop+1
    END if
  END if
end do
do j=jtot,1,-1
  IF(what(j)==w_operator) then
    IF(op(j)==greaterthan.or.op(j)==lessthan) then
      res%operation(jop)%op = op(j)
      ja=j/2
      do WHILE(.not.l_res(ja))
        ja=ja-1
      end do
      res%operation(jop)%a = ja
      res%operation(jop)%b = j/2+1
      l_res(j/2+1)=.false.
      jop = jop+1
    END if
  END if
end do

DEALLOCATE(l_res)

return
end subroutine eval_res_array

subroutine what_is_next(exprin,what,next_l,next_r,op,list_var)
implicit none
CHARACTER(LEN=*), INTENT(IN) :: exprin
INTEGER, INTENT(OUT) :: what,next_l,next_r,op
CHARACTER(*), DIMENSION(:),INTENT(IN), optional :: list_var

CHARACTER(LEN=LEN(exprin)) :: expr

INTEGER :: ln, i, asc, istart,ln2,ierror,jv

ln = LEN_TRIM(ADJUSTL(exprin))
expr = TRIM(ADJUSTL(exprin))
istart = LEN(exprin)-LEN(ADJUSTL(exprin))
next_l = istart

i = 1
asc = IACHAR(expr(1:1))
select case (asc)
 ! parenthesis (
  case (40)
    what = w_parenthesis
    next_r = find_close_bracket(expr(i:ln))
 ! alphabet
  case (65:90,97:122)
    ln2 = ln-i+1
    op = find_intrinsic(expr(i:ln),ln2,ierror)
    IF(ierror==0) then
      what = w_intrinsic
      IF( expr(i+ln2:i+ln2) .ne.'(' ) then
        WRITE(*,*) 'Error in intrinsic ',exprin,': missing parenthesis.'
        stop
      END if
      next_l = next_l+ln2
      next_r = ln2+find_close_bracket(expr(i+ln2:ln))
    ELSEIF(PRESENT(list_var)) then
      what = w_variable
      jv = 0
      op = 0
      do WHILE(jv+1<=SIZE(list_var))
        jv=jv+1
        IF(expr(i:i+ln2-1)==TRIM(list_var(jv)(:))) then
          op=jv
        END if
      end do
      IF(op==0) then
        WRITE(*,*) 'Error: variable ',exprin(i:i+ln2-1),' not allowed.'
        stop
      END if
      next_l = next_l+ln2
      next_r = ln2
    else
        WRITE(*,*) 'Error: intrinsic ',exprin(i:i+ln2-1),' does not exists.'
        stop
    END if
 ! number
  case (48:57)
    what = w_number
    next_r = next_l+find_number_end(expr(i:ln))
 ! operator
  case (42,43,45,47,60,62)
    what = w_operator
    op = find_operator(expr(i:i+1))
    IF(op==power) THEN
      next_r = next_l+2
    else
      next_r = next_l+1
    END if
  case default
end select

return
end subroutine what_is_next

function find_close_bracket(a)
implicit none
INTEGER :: find_close_bracket
CHARACTER(*), INTENT(IN) :: a
INTEGER :: level, i, ln
level = 0
ln=LEN(a)
do i = 1, ln
  select case (a(i:i))
    case ('(')
      level=level+1
    case (')')
      level=level-1
    case default
  end select
  IF(level==0) exit
end do
IF(i>ln) then
  WRITE(*,*) 'Error, no closing bracket in expression ',a
  stop
END if

find_close_bracket = i

return
end function find_close_bracket

function find_intrinsic(a,ln,ierror)
implicit none
INTEGER :: find_intrinsic
CHARACTER(*), INTENT(IN) :: a
CHARACTER(LEN=LEN(a)) :: ac
INTEGER, INTENT(IN OUT) :: ln
INTEGER, INTENT(OUT) :: ierror
INTEGER :: i

ierror=0
do i = 1, ln
  select case (IACHAR(a(i:i)))
    case (65:90,97:122) ! alphabet
    case default
      exit
  end select
end do
IF(i>ln+1) then
  WRITE(*,*) 'Error  in expression ',a
  stop
END if

ln = i-1
ac(1:ln) = a(1:ln)

do i = 1, ln
  IF(IACHAR(a(i:i))<97)   ac(i:i)=ACHAR(IACHAR(a(i:i))+32)
END do

select case (ac(1:ln))
  case ('sqrt')
    find_intrinsic = square_root
  case ('sin')
    find_intrinsic = sine
  case ('cos')
    find_intrinsic = cosine
  case ('tan')
    find_intrinsic = tangent
  case ('exp')
    find_intrinsic = exponential
  case ('log','ln')
    find_intrinsic = logarithm
  case ('log10')
    find_intrinsic = logten
  case ('asin')
    find_intrinsic = arcsin
  case ('acos')
    find_intrinsic = arccos
  case ('atan')
    find_intrinsic = arctan
  case ('sinh')
    find_intrinsic = sinhyp
  case ('cosh')
    find_intrinsic = coshyp
  case ('tanh')
    find_intrinsic = tanhyp
  case default
    ierror=1
end select
return
end function find_intrinsic

function find_operator(a)
implicit none
INTEGER :: find_operator
CHARACTER(*), INTENT(IN) :: a

  select case (IACHAR(a(1:1)))
    case (42) ! '*'
     IF(IACHAR(a(2:2))==42) then
      find_operator = power
     else
      find_operator = multiply
     END if
    case (43) ! '+'
      find_operator = plus
    case (45) ! '-'
      find_operator = minus
    case (47) ! '/'
      find_operator = divide
    case (60) ! '<'
      find_operator = lessthan
    case (62) ! '>'
      find_operator = greaterthan
    case default
      WRITE(*,*) 'Error, operator ',a(1:1),' does not exist'
      stop
  end select

  return
end function find_operator

function find_number_end(a)
implicit none
INTEGER :: find_number_end
CHARACTER(*), INTENT(IN) :: a
INTEGER :: i
LOGICAL :: l_point,l_minus,l_e

l_point=.FALSE.
l_minus=.FALSE.
l_e=.FALSE.

do i = 1, len(a)
  select case (IACHAR(a(i:i)))
    case (48:57)
      l_minus = .true.
    case (46)     ! .
      IF(l_point) exit
      l_point = .true.
      l_minus = .true.
    case (45)     ! -
      IF(l_minus) exit
      l_minus = .true.
    case (69,101) ! e
      IF(l_e) exit
      l_e = .true.
      l_minus = .false.
    case default
      exit
  end select
end do

find_number_end = i-1

return
end function find_number_end

function char2real(a)
implicit none
REAL(amrex_real) :: char2real
CHARACTER(*), INTENT(IN) :: a

READ(a(:),*)  char2real

return
end function char2real

FUNCTION eval(a,b,op)
implicit none
REAL(amrex_real) :: eval
REAL(amrex_real), INTENT(IN) :: a, b
INTEGER, INTENT(IN) :: op

select case (op)
  case (plus)
    eval=a+b
  case (minus)
    eval=a-b
  case (multiply)
    eval=a*b
  case (divide)
    eval=a/b
  case (power)
    eval=a**b
  case (greaterthan)
    eval=0.
    if (a > b) eval=1.
  case (lessthan)
    eval=0.
    if (a < b) eval=1.
  case (square_root)
    eval=SQRT(a)
  case (exponential)
    eval=EXP(a)
  case (logarithm)
    eval=LOG(a)
  case (sine)
    eval=SIN(a)
  case (cosine)
    eval=COS(a)
  case (tangent)
    eval=TAN(a)
  case (logten)
    eval=log10(a)
  case (arcsin)
    eval=asin(a)
  case (arccos)
    eval=acos(a)
  case (arctan)
    eval=atan(a)
  case (sinhyp)
    eval=sinh(a)
  case (coshyp)
    eval=cosh(a)
  case (tanhyp)
    eval=tanh(a)
  case default
    WRITE(*,*) 'Error in eval, wrong operator ',op
    stop
end select

return
END function eval


FUNCTION eval_array(ai,bi,op)
implicit none
REAL(amrex_real), dimension(:), pointer :: eval_array
REAL(amrex_real), INTENT(IN), dimension(:), pointer :: ai, bi
REAL(amrex_real), dimension(:), pointer :: a, b
INTEGER, INTENT(IN) :: op
integer::i1,i2,n

i1 = size(ai)
i2 = size(bi)
n = maxval((/i1,i2/))

select case (op)
  case (plus,minus,multiply,divide,power, greaterthan, lessthan)

  if (i1<n) then
	allocate(a(n))
	a(:) = ai(1)
  else
	a => ai
  end if
  if (i2<n) then
	allocate(b(n))
	b(:) = bi(1)
  else
	b => bi
  end if

  case default
    a => ai
    
end select

allocate(eval_array(n))
eval_array=0.

select case (op)
  case (plus)
    eval_array=a+b
  case (minus)
    eval_array=a-b
  case (multiply)
    eval_array=a*b
  case (divide)
    eval_array=a/b
  case (power)
    eval_array=a**b
  case (square_root)
    eval_array=SQRT(a)
  case (exponential)
    eval_array=EXP(a)
  case (logarithm)
    eval_array=LOG(a)
  case (sine)
    eval_array=SIN(a)
  case (cosine)
    eval_array=COS(a)
  case (tangent)
    eval_array=TAN(a)
  case (logten)
    eval_array=log10(a)
  case (arcsin)
    eval_array=asin(a)
  case (arccos)
    eval_array=acos(a)
  case (arctan)
    eval_array=atan(a)
  case (sinhyp)
    eval_array=sinh(a)
  case (coshyp)
    eval_array=cosh(a)
  case (tanhyp)
    eval_array=tanh(a)
  case default
    WRITE(*,*) 'Error in eval_array, wrong operator ',op
    stop
end select

return
END function eval_array

function C_up2low(cline)
CHARACTER(*) :: cline
CHARACTER(LEN=LEN(cline)) :: C_up2low
INTEGER :: ln, i, j
  ln = LEN(cline)
  do i = 1, ln
    j = IACHAR(cline(i:i))
    IF(j<91.and.j>64) then
     C_up2low(i:i)=ACHAR(j+32)
   else
     C_up2low(i:i) = cline(i:i)
   END if
  END do
end function C_up2low

RECURSIVE FUNCTION C_REPLACE(ST,R1,R2) RESULT(C_RES)
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: ST
CHARACTER(LEN=LEN(st)) :: c_res
CHARACTER(LEN=*), INTENT(IN) :: R1
CHARACTER(LEN=*), INTENT(IN) :: R2

INTEGER :: i,ll

   c_res = ST
   i=INDEX(TRIM(st),r1)
   IF(i==0) return
   ll = LEN(r2)
   IF(ll>0) c_res(i:i+ll-1) = r2
   c_res(i+ll:) = c_replace(c_res(i+LEN(r1):),r1,r2)


RETURN
END FUNCTION C_REPLACE

FUNCTION C_nboccur(ST,R1)
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: ST
integer :: c_nboccur
CHARACTER(LEN=*), INTENT(IN) :: R1

INTEGER :: i,is

  C_nboccur=0
  is=0
10  i = INDEX(st(is+1:),r1)
    IF(i/=0) then
      C_nboccur = C_nboccur+1
      is=is+i+LEN(r1)
      GO TO 10
    END if

RETURN
END FUNCTION C_nboccur

function stringlist(strin) RESULT(res)
CHARACTER(*), INTENT(IN) :: strin
INTEGER :: i1, nwords, i, i2

CHARACTER(LEN(strin)+1) :: str
CHARACTER(80), DIMENSION(:), pointer :: res

str=ADJUSTL(strin)//' '
i1=1
nwords = 0
do WHILE(TRIM(str(i1:))/='')
  i2 = 1
  i2 = INDEX(str(i1:),' ')
  do WHILE(i2==1)
    i1 = i1 + 1
    i2 = INDEX(str(i1:),' ')
  end do
  i1=i1+i2
  nwords = nwords + 1
end do
ALLOCATE(res(nwords))
res(:) = ''
i1=1
do i = 1, nwords
  i2 = INDEX(str(i1:),' ')
  do WHILE(i2==1)
    i1 = i1 + 1
    i2 = INDEX(str(i1:),' ')
  end do
  res(i)(1:i2-1) = str(i1:i1+i2-1)
  i1=i1+i2
end do

return
end function stringlist

SUBROUTINE csv2list(strin, strinlen, lofstr)
CHARACTER(*), INTENT(IN) :: strin
INTEGER, INTENT(IN) :: strinlen
INTEGER :: isep, isep_rel, nwords, i, j
INTEGER, DIMENSION(strinlen) :: word_start, word_end
CHARACTER(strinlen) :: str
CHARACTER(len=:), DIMENSION(:), POINTER, INTENT(inout) :: lofstr(:)
nwords = 1
str = strin
isep_rel = INDEX(str, ',')
isep = isep_rel
word_start(nwords) = 1
DO WHILE(isep_rel > 0)
  word_end(nwords) = isep-1
  nwords = nwords + 1
  word_start(nwords) = isep+1
  str = strin(isep+1:strinlen)
  isep_rel = INDEX(str, ',')
  isep = isep_rel + isep
ENDDO
ALLOCATE( CHARACTER(1) :: lofstr(nwords) )
word_end(nwords) = strinlen
! word_start = word_start(:nwords)
! word_end = word_end(:nwords)
DO i=1, nwords
  lofstr(i)(1:word_end(i)-word_start(i)+1) = strin(word_start(i):word_end(i))
ENDDO
END SUBROUTINE csv2list
end module mod_interpret

! ---------------------
! MODULE parser_wrapper
! Wrappers to use the parser defined in MODULE mod_interpret.
MODULE parser_wrapper
USE iso_c_binding
USE amrex_fort_module, only : amrex_real
USE mod_interpret

IMPLICIT NONE

CONTAINS

! FUNCTION parser_initialize_function RESULT my_index_res
! Initialize a res_type and assign it a unique identifier.
! INPUTS:
!> instr_func   : CHAR* for a mathematical expression
!>                e.g. instr_func = "3*cos(x)+y".
!> instr_var    : CHAR* for a comma-separated list of variables in instr_func
!>                e.g. instr_var = "x,y".
! OUTPUT:
!> my_index_res : INT parser index. Necessary if this module is used for several
!>                parsers (e.g. one for the electron density, one for the ion
!>                density, one for the laser profile etc.).
FUNCTION parser_initialize_function(instr_func, instr_var) RESULT(my_index_res) & 
    bind(c,name='parser_initialize_function')
  CHARACTER(kind=c_char), DIMENSION(*), INTENT(IN) :: instr_func
  CHARACTER(kind=c_char), DIMENSION(*), INTENT(IN) :: instr_var
  INTEGER :: length_func, length_var
  INTEGER :: i
  CHARACTER(kind=c_char, len=:), DIMENSION(:), POINTER :: lofstr(:)
  CHARACTER(kind=c_char, len=120) :: str_func, str_var
  REAL(amrex_real), DIMENSION(1) :: list_var
  INTEGER :: my_index_res

  nb_res = nb_res + 1
  my_index_res = nb_res
  IF (nb_res>10) THEN
    WRITE(*,*) 'Parser error: cannot have more than 10 parsers.'
    STOP
  ENDIF
  
  ! Get length and reformat both inputs before applying the parser.
  length_func=0
  DO
    IF (instr_func(length_func+1) == C_NULL_CHAR) EXIT
    length_func = length_func + 1
  ENDDO
  DO i=1, length_func
    str_func(i:i) = instr_func(i)
  ENDDO
  str_func = str_func(1:length_func)
  length_var=0
  DO
    IF (instr_var(length_var+1) == C_NULL_CHAR) EXIT
    length_var = length_var + 1
  ENDDO
  DO i=1, length_var
    str_var(i:i) = instr_var(i)
  ENDDO
  str_var = str_var(1:length_var)
  ! Convert variable list from csv string to list of strings.
  CALL csv2list(str_var, length_var, lofstr)

  ! Initialize the res object
  CALL eval_res(table_of_res(my_index_res), str_func, list_var=lofstr)

END FUNCTION parser_initialize_function

! FUNCTION parser_evaluate_function RESULT out
! Evaluate parsed function
! INPUTS:
!> list_var     : REAL* array of values for variables
!>                e.g. list_var = (/3.14_8,2._8/).
!> nvar         : INT number of variables
!>                e.g. nvar = 2.
!> my_index_res : INT index of the res_type object in table table_of_res.
! OUTPUT:
!> out          : REAL Result.
FUNCTION parser_evaluate_function(list_var, nvar, my_index_res) result(out) & 
    bind(c,name='parser_evaluate_function')
  INTEGER, VALUE, INTENT(IN) :: nvar, my_index_res
  REAL(amrex_real), INTENT(IN) :: list_var(1:nvar)
  REAL(amrex_real) :: out
  ! Evaluate parsed function in table_of_res(my_index_res), of type res_type
  out = calc_res(table_of_res(my_index_res),list_var)

END FUNCTION parser_evaluate_function

END MODULE parser_wrapper
