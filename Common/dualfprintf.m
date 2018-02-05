function dualfprintf(varagrin)
fid = varagrin{1};
message = varagrin{2} ;
if (nargin == 3)
    val = varagrin{3};
fprintf(fid,message,val);
fprintf(message,val);
else
  fprintf(fid,message);
fprintf(message);  
end