
# fixing colors
sed -i 's+\\{+{+g' $1 
sed -i 's+\\}+}+g' $1 
sed -i 's+color+\\color+g' $1 
# fixing model names
sed -i 's+alpha+\\alpha+g' $1 
sed -i 's+sigma+\\sigma+g' $1 
sed -i 's+tau+\\tau+g' $1 
sed -i 's+emptyset+\\emptyset+g' $1 
# putting pluses
sed -i 's/\\alpha \\sigma/\\alpha + \\sigma/g' $1 
sed -i 's/\\alpha  \\tau/\\alpha + \\tau/g' $1 
sed -i 's/\\sigma \\tau/\\sigma + \\tau/g' $1 
# changing tabular into array
sed -i 's+tabular+array+g' $1 
sed -i 's+\\begin{array}{rllll}+$\\begin{array}{rllll}+g' $1 
sed -i 's+\\end{array}+\\end{array}$+g' $1 
sed -i 's+\\begin{table}++g' $1 
sed -i 's+\[ht\]++g' $1 
sed -i 's+\\end{table}+ +g' $1 
sed -i 's+\\centering+ +g' $1 
# putting rows and col names
sed -i 's+1 \& 2 \& 3 \& 4+\&\\multicolumn{3}{c}{\\text{Type of data}} +g' $1 
sed -i 's+1 \&   \& A+\\multirow{3}{*}{\\rotatebox{90}{Type of model}}\& A \&+g' $1 
sed -i 's+2 \&+\&+g' $1 
sed -i 's+3 \&+\&+g' $1 
sed -i 's+4 \&+\&+g' $1 
# removing bars
sed -i 's+\\hline++g' $1 
# changing slashes 
sed -i 's+SLASH+\\+g' $1 


