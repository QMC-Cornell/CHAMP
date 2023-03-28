!
!       hash_funcs.f -- a library of hash table management routines
!
!                                      by
!
!                              Herbert J. Bernstein
!                                Bernstein + Sons
!                    P.O. Box 177, Bellport, NY 11713-0177, USA
!                    Phone: 1-516-286-1339, Fax: 1-516-286-1999
!                       email: yaya@bernstein-plus-sons.com
!
!       work on these routines done in part at Brookhaven National
!       Laboratory, under contract to the U.S. Department of Energy
!
!-------------------------------------------------------------------------------
!
!       Routines
!
!       hash_init          Initializes a hash table controlled list
!                          call hash_init(data_structure_args)
!
!       hash_find          Searches for a string in a list
!                          call hash_find(name,data_structure_args,ifind)
!
!       hash_store         Inserts as new string in a list
!                          call hash_store(name,data_structure_args,ifind)
!
!       hash_value         Integer function returns index into hash_list
!                          ih = hash_value(name,hash_length)
!
!       The necessary data_structure_args for these routines are
!          name_list   -- an array of character strings
!                         character*(*) name_list(list_length)
!          chain_list  -- chain pointers for searches
!                         integer chain_list(list_length)
!          list_length -- the size of the list arrays
!                         integer list_length
!          num_list    -- number of entries in the list
!                         integer num_list
!          hash_table  -- the initial hashed pointers
!                         integer hash_table
!          hash_length -- the size of the hash table
!                         integer hash_length
!
!
!       The two remaining arguments are
!          name        -- string to search for
!                         character*(*) name
!          ifind       -- return value, 0 for not found (hash_find)
!                         or list full (hash_store), otherwise
!                         the index in name_list of the entry
!
!       The relationship among the arrays used is:
!
!       hash_table is an array (preferably of a modest prime
!       dimension) which starts containing all zeros, which are
!       replaced by pointers to entries in name_list, based
!       values returned by hash_value ranging from 1 to hash_length.
!       Each name is placed in name_list.  A initial zero is placed
!       in the matching entry in chain_list, when the first entry
!       is made.  When a new entry with the same hash_value must be
!       placed a pointer is inserted into chain_list to hook the
!       values together.
!
        subroutine hash_init(name_list,chain_list,list_length,num_list,hash_table,hash_length)
!
!       initialization routine for a hash table controlled list
!          name_list   -- a list of character strings
!          chain_list  -- chain pointers for searches
!          list_length -- the size of the list arrays
!          num_list    -- number of entries in the list
!          hash_table  -- the initial hashed pointers
!          hash_length -- the size of the hash table
!
           character*(*) name_list(list_length)
           integer hash_length,list_length,num_list,i
           integer chain_list(list_length)
           integer hash_table(hash_length)
           num_list=0
           do i = 1,hash_length
           hash_table(i)=0
           enddo
           return
           end
        subroutine hash_find(name,name_list,chain_list,list_length,num_list,hash_table,hash_length,ifind)

!       search routine for a hash table controlled list
!          name        -- string to find
!          name_list   -- a list of character strings
!          chain_list  -- chain pointers for searches
!          list_length -- the size of the list arrays
!          num_list    -- number of entries in the list
!          hash_table  -- the initial hashed pointers
!          hash_length -- the size of the hash table
!          ifind       -- returned index or 0

           character*(*) name
           integer hash_length
           character*(*) name_list(list_length)
           integer chain_list(list_length)
           integer hash_table(hash_length)
           integer hash_value
           integer ifind,list_length,num_list,ih,ip
           ifind=0
           ih=hash_value(name,hash_length)
           ip=hash_table(ih)
 100       if (ip.eq.0) return
           if (name_list(ip).eq.name) then
             ifind=ip
             return
           else
             ip=chain_list(ip)
             go to 100
           endif
           end
        subroutine hash_store(name,name_list,chain_list,list_length,num_list, hash_table,hash_length,ifind)

!       store routine for a hash table controlled list
!          name        -- string to find
!          name_list   -- a list of character strings
!          chain_list  -- chain pointers for searches
!          list_length -- the size of the list arrays
!          num_list    -- number of entries in list
!          hash_table  -- the initial hashed pointers
!          hash_length -- the size of the hash table
!          ifind       -- index of entry or 0 (table full)

           character*(*) name
           character*(*) name_list(list_length)
           integer hash_length
           integer chain_list(list_length)
           integer hash_table(hash_length)
           integer hash_value
           integer ifind,list_length,num_list,ih,ip,iq
           ifind=0
           ih = hash_value(name,hash_length)
           ip=hash_table(ih)
           iq=0
 100       if (ip.eq.0) go to 200
           if (name_list(ip).eq.name) then
             ifind=ip
             return
           else
             iq=ip
             ip=chain_list(ip)
             go to 100
           endif
!JT 200       if (num_list.lt.list_length) then
 200       if (num_list.eq.list_length) then !JT
            write(6,*) 'hash_store: limit list_length=',list_length, ' reached. Increase it!' ! JT
           endif !JT
             num_list=num_list+1
             name_list(num_list)=name
             chain_list(num_list)=0
             if (iq.eq.0) then
               hash_table(ih)=num_list
             else
               chain_list(iq)=num_list
             endif
             ifind=num_list
             return
!JT           else
!JT          ifind = 0
!JT          return
!JT           endif
           end
      integer function hash_value(name,hash_length)

!     function to return a hash value of string name to fit
!     a hash table of length hash_length
      character*(*) name
      integer hash_length,id,ii,i,ic,lenn
      lenn = len(name)
      hash_value=1
      id = 0
      do ii = 1,lenn
        i = 1+lenn-ii
        ic = ichar(name(i:i))
        if (ic.ge.65) then
          hash_value=mod(hash_value*(ic-64),hash_length)+1
          id = id+1
          if (id.gt.3) return
        endif
      enddo
      return
      end
