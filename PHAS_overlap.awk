#! /usr/bin/awk

# Input1 is siRNA result, input2 is miRNA targets result

NR==FNR{
    if(length(substr($1,index($1,"_Chr")+1))==4){
	chrname=substr($1, length($1))
    }else{
	chrname=substr($1, length($1)-1);
    }
    a[chrname"\t"$2"\t"$3] = chrname"\t"$2"\t"$3;
}
NR>FNR{
    chr=$16
    for(i in a){
	split(a[i], tmp, "\t");
#	print tmp[1]"\t"tmp[2]"\t"tmp[3];
	if(tmp[1]==chr && $17=="w"){
	    for(m=tmp[2]-21*5;m<=tmp[2]+21*5;m+=21){
		if(m==$18){
		    print $0"\t"tmp[1]"\t"tmp[2]"\t"tmp[3]"\t5prime"
		}
	    }
	}
	if(tmp[1]==chr && $17=="w"){
	    for(m=tmp[3]+1-21*5;m<=tmp[3]+1+21*5;m+=21){
		if(m==$18){
		    print $0"\t"tmp[1]"\t"tmp[2]"\t"tmp[3]"\t3prime"
		}
	    }
	}
	if(tmp[1]==chr && $17=="c"){
	    for(m=tmp[2]-2+1-21*5;m<=tmp[3]-2+1+21*5;m+=21){
		if(m==$18){
		    print $0"\t"tmp[1]"\t"tmp[2]"\t"tmp[3]"\t5prime"
		}
	    }
	}
	if(tmp[1]==chr && $17=="c"){
	    for(m=tmp[2]-2-21*5;m<=tmp[3]-2+21*5;m+=21){
		if(m==$18){
		    print $0"\t"tmp[1]"\t"tmp[2]"\t"tmp[3]"\t3prime"
		}
	    }
	}
    }
}
