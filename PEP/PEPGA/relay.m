function [] = relay(TextItem,varargin)
% RELAY (TextItem,VariableItem) : This function outputs to both screen and logfile
% 	Essentially, this command duplicates the work of two fprintf commands, and
% 	checks to see that the disk file is being written correctly.
% INPUTS:
%	TextItem		String	Text to print, including formats of variables as per fprintf command
%	varargin		Any		(Optional) The associated variable to print, as per fprintf command

LogFile1 = fopen('Runlog.txt','a+');

    Num_Bytes_Written = 0;
    if nargin > 1	% For printing output that includes a variable
        try
            Num_Bytes_Written = fprintf (LogFile1,TextItem, varargin{:});
            fprintf(1,TextItem,varargin{:});
        end
    else
        try
            Num_Bytes_Written = fprintf (LogFile1,TextItem);
            fprintf(1,TextItem);
        end
    end
fclose(LogFile1);

end %function