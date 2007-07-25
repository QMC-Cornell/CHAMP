      subroutine systemflush (int)
      integer int
#     if defined (UNDERSCORE)
      call flush_ (int)
#     else
      call flush (int)
#     endif
      return
      end

