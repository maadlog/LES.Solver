"use strict";

Object.prototype.orElse = function(value)
{
	if (this) return this; else return value;
}

var relative_error = function(Xminus1,X)
{
	return numeric.norminf(numeric.sub(Xminus1,X));
};

var get_epsilon = function()
{
	return document.getElementById("epsilon").value;
};


function Matrix3x3(firstRow,secondRow,thirdRow, name = "AnonMatrix", independent_column = [0,0,0])
{
	this.a11 = firstRow[0].orElse(0);
	this.a12 = firstRow[1].orElse(0);
	this.a13 = firstRow[2].orElse(0);

	this.a21 = secondRow[0].orElse(0);
	this.a22 = secondRow[1].orElse(0);
	this.a23 = secondRow[2].orElse(0);

	this.a31 = thirdRow[0].orElse(0);
	this.a32 = thirdRow[1].orElse(0);
	this.a33 = thirdRow[2].orElse(0);

	this.name = name;

	this.b1 = independent_column[0].orElse(0);
	this.b2 = independent_column[1].orElse(0);
	this.b3 = independent_column[2].orElse(0);

	this.b_asArray = function()
	{
		return [this.b1,this.b2,this.b3];
	}

	this.toString = function()
	{
		return ""+ this.name +":\n |"+this.a11+" "+this.a12+" "+this.a13+"| \n"
						 		+" |"+this.a21+" "+this.a22+" "+this.a23+"| \n"
								+" |"+this.a31+" "+this.a32+" "+this.a33+"| \n"
	}

	this.asArray = function() {
		return [[this.a11, this.a12, this.a13],
			[this.a21, this.a22, this.a23],
			[this.a31, this.a32, this.a33]]

	}

	this.diagonal = function()
	{
		return [this.a11,this.a22,this.a33];
	}

	this.is_diagonally_dominant = function()
	{
		var dominant_row1 = Math.abs(this.a11) >= (Math.abs(this.a12) + Math.abs(this.a13));
		   
		var dominant_row2 = Math.abs(this.a22) >= (Math.abs(this.a21) + Math.abs(this.a23));

		var dominant_row3 = Math.abs(this.a33) >= (Math.abs(this.a31) + Math.abs(this.a32));

		return dominant_row3 && dominant_row2 && dominant_row1;
	}

	this.get_related_T = function() 
	{
		return new Matrix3x3(
				[0,-(this.a12/this.a11),-(this.a13/this.a11)],
				[-(this.a21/this.a22),0,-(this.a23/this.a22)],
				[-(this.a31/this.a33),-(this.a32/this.a33),0],
				"Resulting T",
				[this.b1/this.a11,this.b2/this.a22,this.b3/this.a33]
			);
	};

	this.infinity_norm = function() 
	{
		numeric.norminf(this.asArray());
	};

	this.spectral_radius = function() 
	{
		try 
		{
			return numeric.eig(this.asArray())
				.lambda.x.reduce(function(x,y)
				{ 
					return Math.max(Math.abs(x),Math.abs(y))
				});
		}
		catch(err) {
    		return 9999; //Eigenvalues don't converge;
		}
	};

	this.result_to_screen = function(res)
	{
		var message = "";
		for (var i = 0 ; i < res.length ; i++) {
			message += " X"+res[i].index+" : ["+res[i].value[0]+","+res[i].value[1]+","+res[i].value[2]+"] \n";
		}
		alert(message);
	}

	this.solve_gauss_seidel = function()
	{
		var X;
		var _x,_y,_z;
		var Xminus1;

		var MatA = this.get_related_T().asArray();
		var MatC = this.get_related_T().b_asArray();

		var keep_going = true;
		var index = 0;
		var res = [];
		
		_x = 0;
		_y = 0;
		_z = 0;

		X = [_x,_y,_z];

		res.push({ "value":X , "index": index });

		Xminus1 = X;
		var ignore = true;
		while(keep_going)
		{
			if (ignore || (relative_error(Xminus1,X) > get_epsilon()))
			{
				_x = (_y * MatA[0][1]) + (_z * MatA[0][2]) + (MatC[0]);
				_y = (_x * MatA[1][0]) + (_z * MatA[1][2]) + (MatC[1]);
				_z = (_x * MatA[2][0]) + (_y * MatA[2][1]) + (MatC[2]);

				Xminus1 = X;
				X = [_x,_y,_z];
				index = index + 1;
				ignore = false;
				res.push({ "value":X , "index": index });
			}
			else
			{
				keep_going = false;
			}
			
		}

		this.result_to_screen(res);

	};

	this.solve_jacobi = function()
	{
		var MatA = this.get_related_T().asArray();
		var MatC = this.get_related_T().b_asArray();

		var keep_going = true;
		var index = 0;
		var res = [];

		var X = [0,0,0];
		
		res.push({ "value":X , "index": index });

		var Xminus1 = X;
		X = MatC;
		index = 1;

		while(keep_going)
		{
			res.push({ "value":X , "index": index });
			if (relative_error(Xminus1,X) > get_epsilon())
			{
				Xminus1 = X;
				X = numeric.add(numeric.dot(MatA,X),MatC);
				index = index + 1;
			}
			else
			{
				keep_going = false;
			}
		}

		this.result_to_screen(res);
	};
}

function Matrix2x2(firstRow,secondRow, name = "AnonMatrix", independent_column = [0,0])
{
	this.a11 = firstRow[0].orElse(0);
	this.a12 = firstRow[1].orElse(0);

	this.a21 = secondRow[0].orElse(0);
	this.a22 = secondRow[1].orElse(0);

	this.name = name;

	this.b1 = independent_column[0].orElse(0);
	this.b2 = independent_column[1].orElse(0);

	this.b_asArray = function()
	{
		return [this.b1,this.b2];
	}

	this.toString = function()
	{
		return ""+ this.name +":\n |"+this.a11+" "+this.a12+"| \n"
						 		+" |"+this.a21+" "+this.a22+"| \n";
	}

	this.asArray = function() {
		return [[this.a11, this.a12],
			[this.a21, this.a22]]

	}

	this.diagonal = function()
	{
		return [this.a11,this.a22];
	}

	this.is_diagonally_dominant = function()
	{
		var dominant_row1 = Math.abs(this.a11) >= Math.abs(this.a12);
		   
		var dominant_row2 = Math.abs(this.a22) >= Math.abs(this.a21);

		return dominant_row2 && dominant_row1;
	}

	this.get_related_T = function() 
	{
		return new Matrix2x2(
				[0,-(this.a12/this.a11)],
				[-(this.a21/this.a22),0],
				"Resulting T",
				[this.b1/this.a11,this.b2/this.a22]
			);
	};

	this.infinity_norm = function() 
	{
		numeric.norminf(this.asArray());
	};
	
	this.spectral_radius = function() 
	{
		try 
		{
			return numeric.eig(this.asArray())
				.lambda.x.reduce(function(x,y)
				{ 
					return Math.max(Math.abs(x),Math.abs(y))
				});
		}
		catch(err) {
    		return 9999; //Eigenvalues don't converge;
		}
	};

	this.result_to_screen = function(res)
	{
		var message = "";
		for (var i = 0 ; i < res.length ; i++) {
			message += " X"+res[i].index+" : ["+res[i].value[0]+","+res[i].value[1]+"] \n";
		}
		alert(message);
	}

	this.solve_gauss_seidel = function()
	{
		var X;
		var _x,_y;
		var Xminus1;

		var MatA = this.get_related_T().asArray();
		var MatC = this.get_related_T().b_asArray();

		var keep_going = true;
		var index = 0;
		var res = [];
		
		_x = 0;
		_y = 0;

		X = [_x,_y];

		res.push({ "value":X , "index": index });

		Xminus1 = X;
		var ignore = true;
		while(keep_going)
		{
			if (ignore || (relative_error(Xminus1,X) > get_epsilon()))
			{
				_x = (_y * MatA[0][1]) + (MatC[0]);
				_y = (_x * MatA[1][0]) + (MatC[1]);

				Xminus1 = X;
				X = [_x,_y];
				index = index + 1;
				ignore = false;
				res.push({ "value":X , "index": index });
			}
			else
			{
				keep_going = false;
			}
			
		}

		this.result_to_screen(res);

	};

	this.solve_jacobi = function()
	{
		var MatA = this.get_related_T().asArray();
		var MatC = this.get_related_T().b_asArray();

		var keep_going = true;
		var index = 0;
		var res = [];

		var X = [0,0];
		
		res.push({ "value":X , "index": index });

		var Xminus1 = X;
		X = MatC;
		index = 1;

		while(keep_going)
		{
			res.push({ "value":X , "index": index });
			if (relative_error(Xminus1,X) > get_epsilon())
			{
				Xminus1 = X;
				X = numeric.add(numeric.dot(MatA,X),MatC);
				index = index + 1;
			}
			else
			{
				keep_going = false;
			}
		}

		this.result_to_screen(res);
	};
}



var _3x3 = true;

var toggle_3x3 = function()
{
	_3x3 = !_3x3;

	if(_3x3)
	{
		document.getElementById("a13").style.display = "inline";
		document.getElementById("a23").style.display = "inline";
		document.getElementById("a31").style.display = "inline";
		document.getElementById("a32").style.display = "inline";
		document.getElementById("a33").style.display = "inline";
		document.getElementById("b3").style.display	 = "inline";
	}
	else
	{
		document.getElementById("a13").style.display = "none";
		document.getElementById("a23").style.display = "none";
		document.getElementById("a31").style.display = "none";
		document.getElementById("a32").style.display = "none";
		document.getElementById("a33").style.display = "none";
		document.getElementById("b3").style.display	 = "none";	
	}
};

var get_a = function(){

	if (_3x3)
	{
		return new Matrix3x3(
			[
				document.getElementById("a11").value,
				document.getElementById("a12").value,
				document.getElementById("a13").value
			],
			[
				document.getElementById("a21").value,
				document.getElementById("a22").value,
				document.getElementById("a23").value
			],
			[
				document.getElementById("a31").value,
				document.getElementById("a32").value,
				document.getElementById("a33").value
			],
			"Matrix A",
			[
				document.getElementById("b1").value,
				document.getElementById("b2").value,
				document.getElementById("b3").value	
			]
			);
	} 
	else
	{
		return new Matrix2x2
		(
			[
				document.getElementById("a11").value,
				document.getElementById("a12").value
			],
			[
				document.getElementById("a21").value,
				document.getElementById("a22").value
			],
			"Matrix A",
			[
				document.getElementById("b1").value,
				document.getElementById("b2").value
			]
		)
	}

}


var	validate_a = function(matrix_a)
{
	document.getElementById("valid_diag_not_zero").style.display = "none";
	document.getElementById("valid_diag_dominant").style.display = "none";
	document.getElementById("valid_norm_of_T").style.display = "none";
	document.getElementById("valid_spectral_radius").style.display = "none";

	if ( matrix_a.diagonal().includes(0) )
	{
		return false;
	}
	else
	{
		document.getElementById("valid_diag_not_zero").style.display = "inline";

		if ( matrix_a.is_diagonally_dominant() )
		{
			document.getElementById("valid_diag_dominant").style.display = "inline";
			return true;
		}
		else 
		{

			if ( matrix_a.get_related_T().infinity_norm() < 1 )
			{
				document.getElementById("valid_norm_of_T").style.display = "inline";
				return true;
			}
			else
			{
				if ( matrix_a.get_related_T().spectral_radius() < 1 )
				{
					document.getElementById("valid_spectral_radius").style.display = "inline";
					return true;
				}
				else
				{
					return false;
				}
			}
		}
	}
};

var solve_a = function()
{
	var matrix_a = get_a();
	if (validate_a(matrix_a))
	{
		
		var jacobi = document.getElementById("jacobi").checked;

		var method;

		if (jacobi)
			method = "Jacobi"
		else
			method = "Gauss-Seidel"

		alert("Solving "+ method +" For: \n" + matrix_a.toString());
		
		if (jacobi)
		{
			matrix_a.solve_jacobi();
		}
		else
		{
			matrix_a.solve_gauss_seidel();
		}
	}
	else
	{
		alert("Matrix didn't meet any convergence criteria");
	}
}



