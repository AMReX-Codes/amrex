module interpolate_module

  contains

      function interpolate(r, npts_model, model_r, model_var)
      
!     given the array of model coordinates (model_r), and variable (model_var),
!     find the value of model_var at point r using linear interpolation.
!     Eventually, we can do something fancier here.
      
      double precision, intent(in   ) :: r
      integer         , intent(in   ) :: npts_model
      double precision, intent(in   ) :: model_r(npts_model), model_var(npts_model)
      double precision                :: interpolate

      ! Local variables
      integer                         :: i, id
      double precision                :: slope,minvar,maxvar
      
!     find the location in the coordinate array where we want to interpolate
      do i = 1, npts_model
         if (model_r(i) .ge. r) exit
      enddo
      
      if (i .eq. npts_model+1) i = i-1
      id = i
      
      if (id .eq. 1) then

         slope = (model_var(id+1) - model_var(id))/(model_r(id+1) - model_r(id))
         interpolate = slope*(r - model_r(id)) + model_var(id)

         ! safety check to make sure interpolate lies within the bounding points
         minvar = min(model_var(id+1),model_var(id))
         maxvar = max(model_var(id+1),model_var(id))
         interpolate = max(interpolate,minvar)
         interpolate = min(interpolate,maxvar)
         
      else if (id .eq. npts_model) then

         slope = (model_var(id) - model_var(id-1))/(model_r(id) - model_r(id-1))
         interpolate = slope*(r - model_r(id)) + model_var(id)

         ! safety check to make sure interpolate lies within the bounding points
         minvar = min(model_var(id),model_var(id-1))
         maxvar = max(model_var(id),model_var(id-1))
         interpolate = max(interpolate,minvar)
         interpolate = min(interpolate,maxvar)

      else

         if (r .ge. model_r(id)) then

            slope = (model_var(id+1) - model_var(id))/(model_r(id+1) - model_r(id))
            interpolate = slope*(r - model_r(id)) + model_var(id)

            ! safety check to make sure interpolate lies within the bounding points
            minvar = min(model_var(id+1),model_var(id))
            maxvar = max(model_var(id+1),model_var(id))
            interpolate = max(interpolate,minvar)
            interpolate = min(interpolate,maxvar)

         else

            slope = (model_var(id) - model_var(id-1))/(model_r(id) - model_r(id-1))
            interpolate = slope*(r - model_r(id)) + model_var(id)

            ! safety check to make sure interpolate lies within the bounding points
            minvar = min(model_var(id),model_var(id-1))
            maxvar = max(model_var(id),model_var(id-1))
            interpolate = max(interpolate,minvar)
            interpolate = min(interpolate,maxvar)

          end if

      endif

      end function interpolate

end module interpolate_module
