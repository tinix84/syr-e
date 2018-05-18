## Copyright (C) 2000-2015 Paul Kienzle
## Copyright (C) 2008-2009 Jaroslav Hajek
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {} {} unique_oct (@var{x})
## @deftypefnx {} {} unique_oct (@var{x}, "rows")
## @deftypefnx {} {[@var{y}, @var{i}, @var{j}] =} unique_oct (@dots{})
## @deftypefnx {} {[@var{y}, @var{i}, @var{j}] =} unique_oct (@dots{}, "first")
## @deftypefnx {} {[@var{y}, @var{i}, @var{j}] =} unique_oct (@dots{}, "last")
## @deftypefnx {} {[@var{y}, @var{i}, @var{j}] =} unique_oct (@dots{}, "stable")
## @deftypefnx {} {[@var{y}, @var{i}, @var{j}] =} unique_oct (@dots{}, "sorted")
## @deftypefnx {} {[@var{y}, @var{i}, @var{j}] =} unique_oct (@dots{}, "stable")
## @deftypefnx {} {[@var{y}, @var{i}, @var{j}] =} unique_oct (@dots{}, "legacy")
## Return the unique_oct elements of @var{x}, sorted in ascending order by default.
##
## If the input @var{x} is a column vector then return a column vector;
## Otherwise, return a row vector.  @var{x} may also be a cell array of strings.
##
## If the optional argument @qcode{"rows"} is given then return the unique_oct
## rows of @var{x}, sorted in ascending order by default.  The input must be a
## 2-D matrix to use this option.
##
## If requested, return index vectors @var{i} and @var{j} such that
## @code{@var{y} = @var{x}(@var{i})} and @code{@var{x} = @var{y}(@var{j})}.
## These are row vectors if @qcode{"legacy"} is specified, or columns otherwise.
##
## The option @qcode{"stable"} returns a vector whose elements appear in the
## same order as they were in @var{x}, whereas @qcode{"sorted"} returns the
## elements in increasing order.  The default is @qcode{"sorted"}.
##
## Additionally, one of @qcode{"first"} or
## @qcode{"last"} may be given as an input.  If @qcode{"last"} is specified,
## return the highest possible indices in @var{i} and for @qcode{"stable"}
## ordering, otherwise, if @qcode{"first"} is specified, return the lowest.
## The default is @qcode{"last"} if @qcode{"legacy"} is set, or @qcode{"first"}
## otherwise.
## @seealso{union, intersect, setdiff, setxor, ismember}
## @end deftypefn

function [y, i, j] = unique_oct (x, varargin)

  if (nargin < 1)
    print_usage ();
  elseif (! (isnumeric (x) || islogical (x) || ischar (x) || iscellstr (x)))
    error ("unique_oct: X must be an array or cell array of strings");
  endif

  if (nargin > 1)
    ## parse options
    if (! iscellstr (varargin))
      error ("unique_oct: options must be strings");
    endif

    optrows   = any (strcmp ("rows", varargin));
    optfirst  = any (strcmp ("first", varargin));
    optlast   = any (strcmp ("last", varargin));
    optstable = any (strcmp ("stable", varargin));
    optsorted = any (strcmp ("sorted", varargin));
    optlegacy = any (strcmp ("legacy", varargin));
    if (optfirst && optlast)
      error ('unique_oct: cannot specify both "first" and "last"');
    elseif (optstable && optsorted)
      error ('unique_oct: cannot specify both "stable" and "sorted"');
    elseif (optfirst + optlast + optrows + optstable + optsorted + optlegacy
            != nargin-1)
      known_options = {"rows", "first", "last", "stable", "sorted", "legacy"};
      for i = 2:nargin
        if ! any (strcmp (varargin{i}, known_options))
          error ("unique_oct: invalid option %s", varargin{i});
        endif
      endfor
    endif

    if (optrows && iscellstr (x))
      warning ('unique_oct: "rows" is ignored for cell arrays');
      optrows = false;
    endif
    
    if (! optlegacy && ! optlast)
      optfirst = true;
    end
  else
    optrows = false;
    optfirst = true;            # new default in Matlab
    optstable = false;
    optlegacy = false;
  endif

  ## FIXME: The operations
  ##
  ##   match = (y(1:n-1) == y(2:n));
  ##   y(idx) = [];
  ##
  ## are very slow on sparse matrices.  Until they are fixed to be as
  ## fast as for full matrices, operate on the nonzero elements of the
  ## sparse array as long as we are not operating on rows.

  if (issparse (x) && ! optrows && nargout <= 1)
    if (nnz (x) < numel (x))
      y = unique_oct ([0; nonzeros(x)], varargin{:});
    else
      ## Corner case where sparse matrix is actually full
      y = unique_oct (full (x), varargin{:});
    endif
    return;
  endif

  if (optrows)
    n = rows (x);
    dim = 1;
  else
    n = numel (x);
    dim = (rows (x) == 1) + 1;
  endif

  y = x;
  ## Special cases 0 and 1
  if (n == 0)
    if (! optrows && isempty (x) && any (size (x)))
      if (iscellstr (y))
        y = cell (0, 1);
      else
        y = zeros (0, 1, class (y));
      endif
    endif
    i = j = [];
    return;
  elseif (n == 1)
    i = j = 1;
    return;
  endif

  if (optrows)
    if (nargout > 1 || optstable)
      [y, ii] = sortrows (y);
    else
      y = sortrows (y);
    endif
    match = all (y(1:n-1,:) == y(2:n,:), 2);
    y(match,:) = [];
  else
    if (! isvector (y))
      y = y(:);
    endif
    if (nargout > 1 || optstable)
      [y, ii] = sort (y);
    else
      y = sort (y);
    endif
    if (iscellstr (y))
      match = strcmp (y(1:n-1), y(2:n));
    else
      match = (y(1:n-1) == y(2:n));
    endif
    y(match) = [];
  endif
  
  ## Calculate output 2 first, since 3 needs perm if optstable
  if (isargout (2) || optstable)
    idx = find (match);
    if (optfirst)
      idx += 1;   # in-place is faster than other forms of increment
    endif
    i = ii;
    if (! isargout (3))
      clear ii            # allow operations on i to avoid copy
    end
    i(idx) = [];
    
    if (optstable)
      [i perm] = sort (i);
      if (optrows)
        y = x(i,:);
      else
        y = x(i);
      endif
    endif
    
    if (! optlegacy)       # New default in Matlab: return column
      i = i(:);
    endif
  endif

  if (isargout (3))
    j = ii;
    positions = cumsum ([1; ! match(:)]);
    
    if (optstable)
      j(ii) = perm(positions);
    else
      j(ii) = positions;
    endif
    
    if (! optlegacy)
      j = j(:);
    endif
  endif
  
endfunction


%!assert (unique_oct ([1 1 2; 1 2 1; 1 1 2]), [1;2])
%!assert (unique_oct ([1 1 2; 1 0 1; 1 1 2],"rows"), [1 0 1; 1 1 2])
%!assert (unique_oct ([]), [])
%!assert (unique_oct ([1]), [1])
%!assert (unique_oct ([1 2]), [1 2])
%!assert (unique_oct ([1;2]), [1;2])
%!assert (unique_oct ([1,NaN,Inf,NaN,Inf]), [1,Inf,NaN,NaN])
%!assert (unique_oct ({"Foo","Bar","Foo"}), {"Bar","Foo"})
%!assert (unique_oct ({"Foo","Bar","FooBar"}'), {"Bar","Foo","FooBar"}')
%!assert (unique_oct (zeros (1,0)), zeros (0,1))
%!assert (unique_oct (zeros (1,0), "rows"), zeros (1,0))
%!assert (unique_oct (cell (1,0)), cell (0,1))
%!assert (unique_oct ({}), {})
%!assert (unique_oct ([1,2,2,3,2,4], "rows"), [1,2,2,3,2,4])
%!assert (unique_oct ([1,2,2,3,2,4]), [1,2,3,4])
%!assert (unique_oct ([1,2,2,3,2,4]', "rows"), [1,2,3,4]')
%!assert (unique_oct (sparse ([2,0;2,0])), [0,2]')
%!assert (unique_oct (sparse ([1,2;2,3])), [1,2,3]')
%!assert (unique_oct ([1,2,2,3,2,4]', "rows"), [1,2,3,4]')
%!assert (unique_oct (single ([1,2,2,3,2,4]), "rows"), single ([1,2,2,3,2,4]))
%!assert (unique_oct (single ([1,2,2,3,2,4])), single ([1,2,3,4]))
%!assert (unique_oct (single ([1,2,2,3,2,4]'), "rows"), single ([1,2,3,4]'))
%!assert (unique_oct (uint8 ([1,2,2,3,2,4]), "rows"), uint8 ([1,2,2,3,2,4]))
%!assert (unique_oct (uint8 ([1,2,2,3,2,4])), uint8 ([1,2,3,4]))
%!assert (unique_oct (uint8 ([1,2,2,3,2,4]'), "rows"), uint8 ([1,2,3,4]'))

%!test
%! [a,i,j] = unique_oct ([1,1,2,3,3,3,4]);
%! assert (a, [1,2,3,4]);
%! assert (i, [2,3,6,7]);
%! assert (j, [1,1,2,3,3,3,4]);
%!
%!test
%! [a,i,j] = unique_oct ([1,1,2,3,3,3,4]', "first");
%! assert (a, [1,2,3,4]');
%! assert (i, [1,3,4,7]');
%! assert (j, [1,1,2,3,3,3,4]');
%!
%!test
%! [a,i,j] = unique_oct ({"z"; "z"; "z"});
%! assert (a, {"z"});
%! assert (i, [3]');
%! assert (j, [1;1;1]);
%!
%!test
%! A = [1,2,3;1,2,3];
%! [a,i,j] = unique_oct (A, "rows");
%! assert (a, [1,2,3]);
%! assert (A(i,:), a);
%! assert (a(j,:), A);

## Test input validation
%!error unique_oct ()
%!error <X must be an array or cell array of strings> unique_oct ({1})
%!error <options must be strings> unique_oct (1, 2)
%!error <cannot specify both "first" and "last"> unique_oct (1, "first", "last")
%!error <invalid option> unique_oct (1, "middle")
%!error <invalid option> unique_oct ({"a", "b", "c"}, "UnknownOption")
%!error <invalid option> unique_oct ({"a", "b", "c"}, "UnknownOption1", "UnknownOption2")
%!error <invalid option> unique_oct ({"a", "b", "c"}, "rows", "UnknownOption2")
%!error <invalid option> unique_oct ({"a", "b", "c"}, "UnknownOption1", "last")
%!warning <"rows" is ignored for cell arrays> unique_oct ({"1"}, "rows");

