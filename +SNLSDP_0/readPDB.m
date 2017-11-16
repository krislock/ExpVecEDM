%%******************************************************************
%% This function is used to extract the coordinate info 
%% from a '.pdb' structure file. Two parameters are 
%% needed for this function: 
%% PDB_file is supposed to be a '.pdb' file.
%% coord_file is supposed to a '.dat' file.
%% The function returns a value in ret_info. 
%% ret_info = 1 if the function suceeds and 0 otherwise.
%%******************************************************************

  function [coords] = readPDB(PDB_file)

%% Count the number of lines of the file1.

  fid = fopen (PDB_file,'r');
  numlines = 0;  ENDT = ' ';
  while ~strcmp(ENDT,'END')
     numlines = numlines + 1;
     A (numlines, 1:79) = fscanf (fid, '%c', [1,79]);
     A (numlines, 80) = fscanf (fid, '%c\n', [1,1]);
     ENDT = A (numlines, 1:3);
  end 
  fclose (fid);

%% Store the coordinates of atoms into coords matrix. row = 0;
   row = 0;
   for i = 1 : numlines
      if  strcmp(A (i,1:4),'ATOM')
         row = row + 1;
         coords (1,row) = sscanf(A(i,31:38),'%f');
         coords (2,row) = sscanf(A(i,39:46),'%f');
         coords (3,row) = sscanf(A(i,47:54),'%f');
      end
   end
%%******************************************************************
