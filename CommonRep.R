cat("\014");
rm(list = objects());

# a

frepcomgen <- function(n, m)
{
	Mat <- matrix(NA, n + 2, m + 2);

	Mat[n + 2, m + 2] <- 1;

	# Puts the values of X on the first column
	for (i in 1 : n)
	{
		Mat[i + 1, 1] <- i;
	}

	# Puts the values of Y on the first line
	for (i in 1 : m)
	{
		Mat[1, i + 1] <- i;
	}

	SampleSize <- 10 ^ 6;

	FrequencyX <- as.numeric(table(sample(1 : n, SampleSize, TRUE)));
	FrequencyY <- as.numeric(table(sample(1 : m, SampleSize, TRUE)));

	# Calculates the probability for the elements of X
	for (i in 1 : n)
	{
		Mat[i + 1, m + 2] <- FrequencyX[i] / SampleSize;
	}

	# Calculates the probability for the elements of Y
	for (i in 1 : m)
	{
		Mat[n + 2, i + 1] <- FrequencyY[i] / SampleSize;
	}

	# Calculates PI(i, j)
	for (i in 2 : (n + 1))
	{
		for (j in 2 : (m + 1))
		{
			Mat[i, j] <- Mat[i, m + 2] * Mat[n + 2, j];
		}
	}

	CountComplete <- n + 1 + m + 1;

	while (CountComplete != 0)
	{
		CompleteLines <- c();
		CompleteColumns <- c();

		# looks for complete lines
		for (i in 2 : (n + 2))
		{
			ExistsNA <- FALSE;

			for (j in 2 : (m + 2))
			{
				if (is.na(Mat[i, j]))
				{
					ExistsNA <- TRUE;
					break;
				}
			}

			if (!ExistsNA)
			{
				CompleteLines[length(CompleteLines) + 1] <- i;
			}
		}

		# looks for complete columns
		for (j in 2 : (m + 2))
		{
			ExistsNA <- FALSE;

			for (i in 2 : (n + 2))
			{
				if (is.na(Mat[i, j]))
				{
					ExistsNA <- TRUE;
					break;
				}
			}

			if (!ExistsNA)
			{
				CompleteColumns[length(CompleteColumns) + 1] <- j;
			}
		}

		CountComplete <- length(CompleteLines) + length(CompleteColumns);

		if (CountComplete == 0)
		{
			next;
		}

		if (length(CompleteLines) == 1)
		{
			CompleteLines[2] <- CompleteLines[1];
		}

		if (length(CompleteColumns) == 1)
		{
			CompleteColumns[2] <- CompleteColumns[1];
		}

		# Chooses a random element from a random complete column and makes it NA
		if (length(CompleteLines) == 0)
		{
			Line <- sample(2 : (n + 2), 1);
			Column <- sample(CompleteColumns, 1);

			Mat[Line, Column] <- NA;
		}
		# Chooses a random element from a random complete line and makes it NA
		else if (length(CompleteColumns) == 0)
		{
			Line<- sample(CompleteLines, 1);
			Column <- sample(2 : (m + 2), 1);

			Mat[Line, Column] <- NA;
		}
		# chooses a random complete line or a random complete column
		# chooses a random element on the random line/column and makes it NA
		else
		{
			if (sample(c(0, 1), 1) == 0)
			{
				Line<- sample(CompleteLines, 1);
				Column <- sample(2 : (m + 2), 1);

				Mat[Line, Column] <- NA;
			}
			else
			{
				Line <- sample(2 : (n + 2), 1);
				Column <- sample(CompleteColumns, 1);
	
				Mat[Line, Column] <- NA;
			}
		}
		# recomputes the number of complete lines and columns
		CompleteLines <- c();
		CompleteColumns <- c();

		for (i in 2 : (n + 2))
		{
			ExistsNA <- FALSE;

			for (j in 2 : (m + 2))
			{
				if (is.na(Mat[i, j]))
				{
					ExistsNA <- TRUE;
					break;
				}
			}

			if (!ExistsNA)
			{
				CompleteLines[length(CompleteLines) + 1] <- i;
			}
		}

		for (j in 2 : (m + 2))
		{
			ExistsNA <- FALSE;

			for (i in 2 : (n + 2))
			{
				if (is.na(Mat[i, j]))
				{
					ExistsNA <- TRUE;
					break;
				}
			}

			if (!ExistsNA)
			{
				CompleteColumns[length(CompleteColumns) + 1] <- j;
			}
		}

		CountComplete <- length(CompleteLines) + length(CompleteColumns);
	}

	return (Mat);
};

Mat <- frepcomgen(3, 4);
print(Mat);

# b

fcomplrepcom <- function(Mat)
{
	n <- nrow(Mat) - 2;
	m <- ncol(Mat) - 2;

	# counts holes in the matrix
	CountNA <- 0;
	for (i in 2 : (n + 2))
	{
		for (j in 2 : (m + 2))
		{
			if (is.na(Mat[i, j]))
			{
				CountNA <- CountNA + 1;
			}
		}
	}

	if (is.na(Mat[n + 2, m + 2]))
	{
		CountNA <- CountNA - 1;
		Mat[n + 2, m + 2] <- 1;
	}

	while (CountNA != 0)
	{
		Found <- FALSE;

		for (i in 2 : (n + 2))
		{
			AuxCount <- 0;
			AuxJ <- NA;
			# counts the number of holes and the current line
			for (j in 2 : (m + 2))
			{
				if (is.na(Mat[i, j]))
				{
					AuxCount <- AuxCount + 1;
					AuxJ <- j;
				}
			}

			# when there's only 1 number left on the current line it completes the line
			if (AuxCount == 1)
			{
				# when the hole is the last element on the line
				# the probability is the sum of the other probabilities on the line
				if (AuxJ == m + 2)
				{
					Sum <- 0;

					for (j in 2 : (m + 1))
					{
						Sum <- Sum + Mat[i, j];
					}

					Mat[i, AuxJ] <- Sum;
				}
				# otherwise the probability of the current element is 
				# the difference between the last element of the line and the other elements
				else
				{
					Sum <- Mat[i, m + 2];

					for (j in 2 : (m + 1))
					{
						if (!is.na(Mat[i, j]))
						{
							Sum <- Sum - Mat[i, j];
						}
					}

					Mat[i, AuxJ] <- Sum;
				}

				CountNA <- CountNA - 1;
				Found <- TRUE;
				break;
			}
		}

		if (Found == TRUE)
		{
			next;
		}

		# Same thing but looks for the columns where there's one hole
		for (j in 2 : (m + 2))
		{
			AuxCount <- 0;
			AuxI <- NA;

			for (i in 2 : (n + 2))
			{
				if (is.na(Mat[i, j]))
				{
					AuxCount <- AuxCount + 1;
					AuxI <- i;
				}
			}

			if (AuxCount == 1)
			{
				if (AuxI == n + 2)
				{
					Sum <- 0;

					for (i in 2 : (n + 1))
					{
						Sum <- Sum + Mat[i, j];
					}

					Mat[AuxI, j] <- Sum;
				}
				else
				{
					Sum <- Mat[n + 2, j];

					for (i in 2 : (n + 1))
					{
						if (!is.na(Mat[i, j]))
						{
							Sum <- Sum - Mat[i, j];
						}
					}

					Mat[AuxI, j] <- Sum;
				}

				CountNA <- CountNA - 1;
				Found <- TRUE;
				break;
			}
		}

		if (Found == TRUE)
		{
			next;
		}

		# If there's no column nor line with only one element we set the probability
		# of one hole with 0 and we startover.
		for (i in 2 : (n + 1))
		{
			for (j in 2 : (m + 1))
			{
				if (is.na(Mat[i, j]))
				{
					Mat[i, j] <- 0;
					Found <- TRUE;
					break;
				}
			}

			if (Found == TRUE)
			{
				break;
			}
		}

		CountNA <- CountNA - 1;
	}

	return (Mat);
};

Mat <- fcomplrepcom(Mat);
print(Mat);

# c
frepmarginal <- function(Mat)
{
	n <- nrow(Mat) - 2;
	m <- ncol(Mat) - 2;

	X <- matrix(NA, 2, n);
	Y <- matrix(NA, 2, m);

	# the first line of X is the first column of the matrix
	# the second line of X is the last column of the matrix
	for (i in 1 : n)
	{
		X[1, i] <- Mat[i + 1, 1];
		X[2, i] <- Mat[i + 1, m + 2];
	}

	# the first line of Y is the first column of the matrix
	# the second line of Y is the last column of the matrix
	for (i in 1 : m)
	{
		Y[1, i] <- Mat[1, i + 1];
		Y[2, i] <- Mat[n + 2, i + 1];
	}

	return (list(X, Y));
};

List <- frepmarginal(Mat);
X <- List[[1]];
Y <- List[[2]];
rm("List");
print(X);
print(Y);

# d
# The expected value of a variable is the sum of
# the product between the current value and the current probability for all values of var.
expected <- function(Var)
{
	Result <- 0;

	for (i in 1 : ncol(Var))
	{
		Result <- Result + Var[1, i] * Var[2, i];
	}

	return (Result);
};

fpropcov <- function(a, b, c, d, Mat)
{
	List <- frepmarginal(Mat);

	X <- List[[1]];
	Y <- List[[2]];

	XX <- matrix(NA, 2, ncol(X));
	XY <- matrix(NA, 2, ncol(X) * ncol(Y));
	YY <- matrix(NA, 2, ncol(Y));

	# X^2 is the values of x squared with the probabilities of x
	for (i in 1 : ncol(X))
	{
		XX[1, i] <- X[1, i] * X[1, i];
		XX[2, i] <- X[2, i];
	}

	# We take all the combinations of values betweeen x and y
	for (i in 1 : ncol(X))
	{
		for (j in 1 : ncol(Y))
		{
			XY[1, i + (j - 1) * ncol(X)] <- Mat[i + 1, 1] * Mat[1, j + 1];
			XY[2, i + (j - 1) * ncol(X)] <- Mat[i + 1, j + 1];
		}
	}

	# Y^2 is the values of y squared with the probabilities of y
	for (i in 1 : ncol(Y))
	{
		YY[1, i] <- Y[1, i] * Y[1, i];
		YY[2, i] <- Y[2, i];
	}

	ExpX <- expected(X);
	ExpY <- expected(Y);
	ExpXX <- expected(XX);
	ExpXY <- expected(XY);
	ExpYY <- expected(YY);

	# we rewrote the formula of covariance as ac * Var(X) + (ad + bc) * cov(X, Y) + bd * Var(Y)
	return (a * c * (ExpXX - ExpX * ExpX) + (a * d + b * c) * (ExpXY - ExpX * ExpY) + b * d * (ExpYY - ExpY * ExpY));
};

print(fpropcov(1, 2, 3, 4, Mat));

# e

fPcond <- function(Mat, XVal, YVal, XDependsOnY)
{
	XIndex <- 1;

	# searches the index of the element XVal in X
	while (Mat[XIndex + 1, 1] != XVal)
	{
		XIndex <- XIndex + 1;
	}

	YIndex <- 1;

	# searches the index of the element YVal in Y
	while (Mat[1, YIndex + 1] != YVal)
	{
		YIndex <- YIndex + 1;
	}

	# P(X = x | Y = y) = P(X = x, Y = y) / P(Y = y)
	if (XDependsOnY == TRUE)
	{
		return (Mat[XIndex + 1, YIndex + 1] / Mat[nrow(Mat), YIndex + 1]);
	}

	
	# P(Y = y | X = x) = P(X = x, Y = y) / P(X = x)
	return (Mat[XIndex + 1, YIndex + 1] / Mat[XIndex + 1, ncol(Mat)]);
};

print(fPcond(Mat, 1, 2, TRUE));
print(fPcond(Mat, 1, 2, FALSE));

# f

fPcomun <- function(Mat, XVal, YVal)
{
	XIndex <- 1;

	# searches the index of the element XVal in X
	while (Mat[XIndex + 1, 1] != XVal)
	{
		XIndex <- XIndex + 1;
	}

	YIndex <- 1;

	# searches the index of the element YVal in Y
	while (Mat[1, YIndex + 1] != YVal)
	{
		YIndex <- YIndex + 1;
	}

	return (Mat[XIndex + 1, YIndex + 1]);
};

print(fPcomun(Mat, 1, 2));

# g

# 1
# Transforms the values of X in the matrix as g(x) = x * 5 + 9
Mat[2 : (nrow(Mat) - 1), 1] <- Mat[2 : (nrow(Mat) - 1), 1] * 5 + 9;
# Transforms the values of Y in the matrix as h(y) = y * (-3) - 9
Mat[1, 2 : (ncol(Mat) - 1)] <- Mat[1, 2 : (ncol(Mat) - 1)] * (-3) - 2;

# Prints cov(X, Y)
print(fpropcov(1, 0, 0, 1, Mat));

# Removes the transformations from the matrix
Mat[2 : (nrow(Mat) - 1), 1] <- (Mat[2 : (nrow(Mat) - 1), 1] - 9) / 5;
Mat[1, 2 : (ncol(Mat) - 1)] <- (Mat[1, 2 : (ncol(Mat) - 1)] + 2) / (-3);

# 2
PX <- 0;
PY <- 0;

# Computes P(Y > 0.3)
for (i in 2 : (ncol(Mat) - 1))
{
	if (0.3 < Mat[1, i])
	{
		PY <- PY + Mat[nrow(Mat), i];
	}
}

# Computes P(0 < X < 0.8, Y < 1.7)
for (i in 2 : (nrow(Mat) - 1))
{
	if (0 >= Mat[i, 1] || Mat[i, 1] >= 0.8)
	{
		next;
	}

	for (j in 2 : (ncol(Mat) - 1))
	{
		if (0.3 >= Mat[1, j])
		{
			next;
		}

		PX <- PX + Mat[i, j];
	}
}


print(PX / PY);

if (exists("i"))
{
	rm("i");
}
if (exists("j"))
{
	rm("j");
}
rm("PX");
rm("PY");

# 3

P <- 0;
# Sums the probabilities where X > 0.2 and Y < 1.7
for (i in 2 : (nrow(Mat) - 1))
{
	if (0.2 >= Mat[i, 1])
	{
		next;
	}

	for (j in 2 : (ncol(Mat) - 1))
	{
		if (Mat[1, j] >= 1.7)
		{
			next;
		}

		P <- P + Mat[i, j];
	}
}

print(P);

if (exists("i"))
{
	rm("i");
}
if (exists("j"))
{
	rm("j");
}
rm("P");

# h

# 1

# Checks if PI(i, j) = p(i) * q(j) for all i, j
fverind <- function(Mat)
{
	for (i in 2 : (nrow(Mat) - 1))
	{
		for (j in 2 : (ncol(Mat) - 1))
		{
			if (abs(Mat[i, j] - Mat[i, ncol(Mat)] * Mat[nrow(Mat), j]) > 1e-8)
			{
				return (FALSE);
			}
		}
	}

	return (TRUE);
};

print(fverind(Mat));

# 2

# if cov(X, Y) is 0 then the variables are not correlated
fvernecor <- function(Mat)
{
	return (abs(fpropcov(1, 0, 0, 1, Mat)) < 1e-8);
};

print(fvernecor(Mat));

# i
# Generates 3D matrix
generate3D <- function(X, Y, Z)
{
	XLen <- ncol(X);
	YLen <- ncol(Y);
	ZLen <- ncol(Z);

	Mat <- array(NA, c(XLen + 2, YLen + 2, ZLen + 2));

	for (i in 1 : XLen)
	{
		# Places values of X on the edge of the cube
		Mat[i + 1, 1, 1] <- X[1, i];
		# Places probability of X = x on the edges of the cube that are paralel to the values of X
		Mat[i + 1, YLen + 2, 1] <- X[2, i];
		Mat[i + 1, 1, ZLen + 2] <- X[2, i];
		Mat[i + 1, YLen + 2, ZLen + 2] <- X[2, i];
	}

	for (i in 1 : YLen)
	{
		# Places values of Y on the edge of the cube
		Mat[1, i + 1, 1] <- Y[1, i];
		# Places probability of Y = y on the edges of the cube that are paralel to the values of Y
		Mat[1, i + 1, ZLen + 2] <- Y[2, i];
		Mat[XLen + 2, i + 1, 1] <- Y[2, i];
		Mat[XLen + 2, i + 1, ZLen + 2] <- Y[2, i];
	}

	for (i in 1 : ZLen)
	{
		# Places values of Z on the edge of the cube
		Mat[1, 1, i + 1] <- Z[1, i];
		# Places probability of Z = z on the edges of the cube that are paralel to the values of Z
		Mat[XLen + 2, 1, i + 1] <- Z[2, i];
		Mat[1, YLen + 2, i + 1] <- Z[2, i];
		Mat[XLen + 2, YLen + 2, i + 1] <- Z[2, i];
	}

	# Places 1 in the oposit corner of the cube

	Mat[XLen + 2, YLen + 2, ZLen + 2] <- 1;

	# Generates the common repartition for XYZ

	for (i in 2 : (XLen + 1))
	{
		for (j in 2 : (YLen + 1))
		{
			for (k in 2 : (ZLen + 1))
			{
				Mat[i, j, k] <- X[2, i - 1] * Y[2, j - 1] * Z[2, k - 1];
			}
		}
	}

	# Generates the common repartition for XY

	for (i in 2 : (XLen + 1))
	{
		for (j in 2 : (YLen + 1))
		{
			Mat[i, j, 1] <- X[2, i - 1] * Y[2, j - 1];
			Mat[i, j, ZLen + 2] <- X[2, i - 1] * Y[2, j - 1];
		}
	}

	# Generates the common repartition for YZ

	for (i in 2 : (YLen + 1))
	{
		for (j in 2 : (ZLen + 1))
		{
			Mat[1, i, j] <- Y[2, i - 1] * Z[2, j - 1];
			Mat[XLen + 2, i, j] <- Y[2, i - 1] * Z[2, j - 1];
		}
	}

	# Generates the common repartition for XZ

	for (i in 2 : (XLen + 1))
	{
		for (j in 2 : (ZLen + 1))
		{
			Mat[i, 1, j] <- X[2, i - 1] * Z[2, j - 1];
			Mat[i, YLen + 2, j] <- X[2, i - 1] * Z[2, j - 1];
		}
	}

	return (Mat);
};

# Generates a 3rd random variable in order to use it in the generate3D function

Z <- matrix(NA, 2, 3);

for (i in 1 : 3)
{
	Z[1, i] <- i;
	Z[2, i] <- 1 / 3;
}

print(generate3D(X, Y, Z));
