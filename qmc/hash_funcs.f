C
C       hash_funcs.f -- a library of hash table management routines
C
C                                      by
C
C                              Herbert J. Bernstein
C                                Bernstein + Sons
C                    P.O. Box 177, Bellport, NY 11713-0177, USA
C                    Phone: 1-516-286-1339, Fax: 1-516-286-1999
C                       email: yaya@bernstein-plus-sons.com
C
C       work on these routines done in part at Brookhaven National
C       Laboratory, under contract to the U.S. Department of Energy
C
C-------------------------------------------------------------------------------
C
C       Routines
C
C       hash_init          Initializes a hash table controlled list
C                          call hash_init(data_structure_args)
C
C       hash_find          Searches for a string in a list
C                          call hash_find(name,data_structure_args,ifind)
C
C       hash_store         Inserts as new string in a list
C                          call hash_store(name,data_structure_args,ifind)
C
C       hash_value         Integer function returns index into hash_list
C                          ih = hash_value(name,hash_length)
C
C       The necessary data_structure_args for these routines are
C          name_list   -- an array of character strings
C                         character*(*) name_list(list_length)
C          chain_list  -- chain pointers for searches
C                         integer chain_list(list_length)
C          list_length -- the size of the list arrays
C                         integer list_length
C          num_list    -- number of entries in the list
C                         integer num_list
C          hash_table  -- the initial hashed pointers
C                         integer hash_table
C          hash_length -- the size of the hash table
C                         integer hash_length
C
C
C       The two remaining arguments are
C          name        -- string to search for
C                         character*(*) name
C          ifind       -- return value, 0 for not found (hash_find)
C                         or list full (hash_store), otherwise
C                         the index in name_list of the entry
C
C       The relationship among the arrays used is:
C
C       hash_table is an array (preferably of a modest prime
C       dimension) which starts containing all zeros, which are
C       replaced by pointers to entries in name_list, based
C       values returned by hash_value ranging from 1 to hash_length.
C       Each name is placed in name_list.  A initial zero is placed
C       in the matching entry in chain_list, when the first entry
C       is made.  When a new entry with the same hash_value must be
C       placed a pointer is inserted into chain_list to hook the
C       values together.
C
        subroutine hash_init(name_list,chain_list,list_length,num_list,
     *                       hash_table,hash_length)
C
C       initialization routine for a hash table controlled list
C          name_list   -- a list of character strings
C          chain_list  -- chain pointers for searches
C          list_length -- the size of the list arrays
C          num_list    -- number of entries in the list
C          hash_table  -- the initial hashed pointers
C          hash_length -- the size of the hash table
C
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
        subroutine
     *  hash_find(name,name_list,chain_list,list_length,num_list,
     *                       hash_table,hash_length,ifind)
C
C       search routine for a hash table controlled list
C          name        -- string to find
C          name_list   -- a list of character strings
C          chain_list  -- chain pointers for searches
C          list_length -- the size of the list arrays
C          num_list    -- number of entries in the list
C          hash_table  -- the initial hashed pointers
C          hash_length -- the size of the hash table
C          ifind       -- returned index or 0
C
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
        subroutine
     *  hash_store(name,name_list,chain_list,list_length,num_list,
     *                       hash_table,hash_length,ifind)
C
C       store routine for a hash table controlled list
C          name        -- string to find
C          name_list   -- a list of character strings
C          chain_list  -- chain pointers for searches
C          list_length -- the size of the list arrays
C          num_list    -- number of entries in list
C          hash_table  -- the initial hashed pointers
C          hash_length -- the size of the hash table
C          ifind       -- index of entry or 0 (table full)
C
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
            write(6,*) 'hash_store: limit list_length=',list_length, !§JT
     *     ' reached. Increase it!'                                  ! JT
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
C
C     function to return a hash value of string name to fit
C     a hash table of length hash_length
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
        
