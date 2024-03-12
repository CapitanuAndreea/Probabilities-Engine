cat("\014");
rm(list = objects());
library(shiny);

# User defined common repartition distribuion

f <- function(x, y)
{
	# Gaussian distribution in 2D

	return (exp(-x ^ 2 - y ^ 2) / 3.1415926535897);
};

# Default transformation for the random variable

Identity <- function(x)
{
	return (x);
};

# Default transformation for the 2D random variable

Identity2 <- function(x, y)
{
	return (x * y);
};

# a

Fubini <- function(f)
{
	# The intervals for the integral

	a <- -10 ^ 1;
	b <- 10 ^ 1;
	c <- -10 ^ 1;
	d <- 10 ^ 1;
	Step <- 0.1;

	# Integral results

	Sum1 <- 0;
	Sum2 <- 0;

	a_copy <- a;
	b_copy <- b;
	c_copy <- c;
	d_copy <- d;

	# Integrating over x

	while (a_copy < b_copy)
	{
		c_copy <- c;

		# Integrating over y

		while (c_copy < d_copy)
		{
			Sum1 <- Sum1 + f(a_copy, c_copy) * Step * Step;

			c_copy <- c_copy + Step;
		}

		a_copy <- a_copy + Step;
	}

	a_copy <- a;
	b_copy <- b;
	c_copy <- c;
	d_copy <- d;

	# Integrating over y

	while (c_copy < d_copy)
	{
		a_copy <- a;

		# Integrating over x

		while (a_copy < b_copy)
		{
			Sum2 <- Sum2 + f(a_copy, c_copy) * Step * Step;

			a_copy <- a_copy + Step;
		}

		c_copy <- c_copy + Step;
	}

	# Checking if both results are equal

	if (abs(Sum1 - Sum2) < 1e-8)
	{
		print(Sum1);
	}
	else
	{
		print("Nu se aplica Fubini");
	}

	return (abs(Sum1 - Sum2) < 1e-8);
};

Fubini(f);

# b

# Defining the GUI layout

GUI <- fluidPage(
	sidebarLayout(
		sidebarPanel(
			sliderInput("XMin", "X Min", min = -10, max = 10, value = -5),
			sliderInput("XMax", "X Max", min = -10, max = 10, value = 5),
			sliderInput("YMin", "Y Min", min = -10, max = 10, value = -5),
			sliderInput("YMax", "Y Max", min = -10, max = 10, value = 5),
			sliderInput("Theta", "Display Theta", min = 0, max = 360, value = 45),
			sliderInput("Phi", "Display Phi", min = 0, max = 360, value = 45)
		),
		mainPanel(
			plotOutput("plot")
		)
	)
);

# Defining the behaviour of the server

Server <- function(input, output)
{
	output$plot <- renderPlot(
		{
			x <- seq(input$XMin, input$XMax, length.out = 30);
			y <- seq(input$YMin, input$YMax, length.out = 30);
			z <- outer(x, y, f);

			persp(x, y, z, theta = input$Theta, phi = input$Phi, col = "green", border = "purple", shade = 1);
		}
	)
};

# Running the app

shinyApp(GUI, Server);

# c

CheckDensity <- function(f)
{
	# The intervals for the integral

	a <- -10 ^ 1;
	b <- 10 ^ 1;
	c <- -10 ^ 1;
	d <- 10 ^ 1;
	Step <- 0.1;

	# The result of the integral

	Sum <- 0;

	# Integrating over x

	while (a < b)
	{
		c_copy <- c;

		# Integrating over y

		while (c_copy < d)
		{
			fVal <- f(a, c_copy);

			# Checking for negative values in the distribution function

			if (fVal < 0)
			{
				return (FALSE);
			}

			Sum <- Sum + fVal * Step * Step;

			c_copy <- c_copy + Step;
		}

		a <- a + Step;
	}

	# Checking if the total probability is 1, as it should be

	return (abs(Sum - 1) < 1e-8);
};

print(CheckDensity(f));

# d

# Saving the random variable in an object. This so called object is a function pointer

XY <- f;

# e

X <- function(x)
{
	# The result of the integral

	Sum <- 0;

	# The intervals for the integral

	a <- -10 ^ 1;
	b <- 10 ^ 1;
	Step <- 0.1;

	# Integrating over y

	while (a < b)
	{
		Sum <- Sum + XY(x, a) * Step;

		a <- a + Step;
	}

	return (Sum);
};

Y <- function(y)
{
	# The result of the integral

	Sum <- 0;

	# The intervals for the integral

	a <- -10 ^ 1;
	b <- 10 ^ 1;
	Step <- 0.1;

	# Integrating over x

	while (a < b)
	{
		Sum <- Sum + XY(a, y) * Step;

		a <- a + Step;
	}

	return (Sum);
};

XCondY <- function(x, y)
{
	# The common density devided by the marginal density of y

	return (XY(x, y) / Y(y));
};

YCondX <- function(y, x)
{
	# The common density devided by the marginal density of x

	return (XY(x, y) / X(x));
};

# f

XRep <- function(x)
{
	# The result of the integral

	Sum <- 0;

	# The intervals of the integral

	a1 <- -10 ^ 1;
	b1 <- x;
	Step <- 0.1;

	# Integrating

	while (a1 < b1)
	{
		Sum <- Sum + X(a1) * Step;

		a1 <- a1 + Step;
	}

	return (Sum);
};

# Abstract multidimensional function used by plot

XRepAbstr <- function(x_list)
{
	Result <- x_list;

	for (i in 1 : length(x_list))
	{
		Result[i] <- XRep(x_list[i]);
	}

	return (Result);
};

YRep <- function(y)
{
	# The result of the integral

	Sum <- 0;

	# The intervals of the integral

	a1 <- -10 ^ 1;
	b1 <- y;
	Step <- 0.1;

	# Integrating

	while (a1 < b1)
	{
		Sum <- Sum + Y(a1) * Step;

		a1 <- a1 + Step;
	}

	return (Sum);
};

# Abstract multidimensional function used by plot

YRepAbstr <- function(y_list)
{
	Result <- y_list;

	for (i in 1 : length(y_list))
	{
		Result[i] <- YRep(y_list[i]);
	}

	return (Result);
};

XYRep <- function(x, y)
{
	# The intervals of the integral

	a <- -10 ^ 1;
	b <- x;
	c <- -10 ^ 1;
	d <- y;
	Step <- 0.1;

	# The result of the integral

	Sum <- 0;

	# Integrating over x

	while (a < b)
	{
		c_copy <- c;

		# Integrating over y

		while (c_copy < d)
		{
			Sum <- Sum + XY(a, c_copy) * Step * Step;

			c_copy <- c_copy + Step;
		}

		a <- a + Step;
	}

	return (Sum);
};

# Abstract multidimensional function used by persp

XYRepAbstr <- function(x_list, y_list)
{
	Result <- x_list;

	for (i in 1 : length(x_list))
	{
		Result[i] <- XYRep(x_list[i], y_list[i]);
	}

	return (Result);
};

# Defining the GUI layout

GUI <- fluidPage(
	sidebarLayout(
		sidebarPanel(
			sliderInput("XMinDens", "X Min Density", min = -10, max = 10, value = -5),
			sliderInput("XMaxDens", "X Max Density", min = -10, max = 10, value = 5),

			sliderInput("XMinRept", "X Min Repartition", min = -10, max = 10, value = -5),
			sliderInput("XMaxRept", "X Max Repartition", min = -10, max = 10, value = 5),

			sliderInput("XY_XMinDens", "XY: X Min Density", min = -10, max = 10, value = -5),
			sliderInput("XY_XMaxDens", "XY: X Max Density", min = -10, max = 10, value = 5),
			sliderInput("XY_YMinDens", "XY: Y Min Density", min = -10, max = 10, value = -5),
			sliderInput("XY_YMaxDens", "XY: Y Max Density", min = -10, max = 10, value = 5),

			sliderInput("XY_XMinRept", "XY: X Min Repartition", min = -10, max = 10, value = -5),
			sliderInput("XY_XMaxRept", "XY: X Max Repartition", min = -10, max = 10, value = 5),
			sliderInput("XY_YMinRept", "XY: Y Min Repartition", min = -10, max = 10, value = -5),
			sliderInput("XY_YMaxRept", "XY: Y Max Repartition", min = -10, max = 10, value = 5),

			sliderInput("Theta", "Display Theta", min = 0, max = 360, value = 45),
			sliderInput("Phi", "Display Phi", min = 0, max = 360, value = 45)
		),
		mainPanel(
			plotOutput("DensityX"),
			plotOutput("RepartitionX"),
			plotOutput("DensityXY"),
			plotOutput("RepartitionXY")
		)
	)
);

# Defining the behaviour of the server

Server <- function(input, output)
{
	output$DensityX <- renderPlot(
		{
			# Density of X
			plot(X, xlim = c(input$XMinDens, input$XMaxDens));
		}
	)
	output$RepartitionX <- renderPlot(
		{
			# Repartition of X
			plot(XRepAbstr, xlim = c(input$XMinRept, input$XMaxRept));
		}
	)
	output$DensityXY <- renderPlot(
		{
			# Density of XY

			x <- seq(input$XY_XMinDens, input$XY_XMaxDens, length.out = 30);
			y <- seq(input$XY_YMinDens, input$XY_YMaxDens, length.out = 30);
			z <- outer(x, y, XY);

			persp(x, y, z, theta = input$Theta, phi = input$Phi, col = "green", border = "purple", shade = 1);
		}
	)
	output$RepartitionXY <- renderPlot(
		{
			# Repartition of XY

			x <- seq(input$XY_XMinRept, input$XY_XMaxRept, length.out = 30);
			y <- seq(input$XY_YMinRept, input$XY_YMaxRept, length.out = 30);
			z <- outer(x, y, XYRepAbstr);

			persp(x, y, z, theta = input$Theta, phi = input$Phi, col = "green", border = "purple", shade = 1);
		}
	)
};

# Running the app

shinyApp(GUI, Server);

# g

Medie <- function(f, power, Transform = Identity)
{
	# The result of the integral

	Sum <- 0;

	# The intervals for the integral

	a <- -10 ^ 1;
	b <- 10 ^ 1;
	Step <- 0.1;

	# Integrating

	while (a < b)
	{
		Sum <- Sum + Transform(a) ^ power * f(a) * Step;

		a <- a + Step;
	}

	return (Sum);
};

Medie2 <- function(f, power, Transform = Identity2)
{
	# The result of the integral

	Sum <- 0;

	# The intervals for the integral

	a <- -10 ^ 1;
	b <- 10 ^ 1;
	c <- -10 ^ 1;
	d <- 10 ^ 1;
	Step <- 0.1;

	# Integrating over x

	while (a < b)
	{
		c_copy <- c;

		# Integrating over y

		while (c_copy < d)
		{
			Sum <- Sum + Transform(a, c_copy) ^ power * f(a, c_copy) * Step * Step;

			c_copy <- c_copy + Step;
		}

		a <- a + Step;
	}

	return (Sum);
};

MomCentr <- function(f, power)
{
	# The result of the integral

	Sum <- 0;

	# The intervals for the integral

	a <- -10 ^ 1;
	b <- 10 ^ 1;
	Step <- 0.1;

	# Computing the expected value for the variable

	Expected <- Medie(f, 1);

	# Integrating

	while (a < b)
	{
		Sum <- Sum + (a - Expected) ^ power * f(a) * Step;

		a <- a + Step;
	}

	return (Sum);
};

MomCentr2 <- function(f, power)
{
	# The result of the integral

	Sum <- 0;

	# The intervals for the integral

	a <- -10 ^ 1;
	b <- 10 ^ 1;
	c <- -10 ^ 1;
	d <- 10 ^ 1;
	Step <- 0.1;

	# Computing the expected value for the product of X and Y using the common repartition density

	Expected <- Medie2(f, 1);

	# Integrating over x

	while (a < b)
	{
		c_copy <- c;

		# Integrating over y

		while (c_copy < d)
		{
			Sum <- Sum + (a * c_copy - Expected) ^ power * f(a, c_copy) * Step * Step;

			c_copy <- c_copy + Step;
		}

		a <- a + Step;
	}

	return (Sum);
};

# Computing the required values

EX <- Medie(X, 1);
print(EX);
EY <- Medie(Y, 1);
print(EY);
EXY <- Medie2(XY, 1);
print(EXY);

EXX <- Medie(X, 2);
EYY <- Medie(Y, 2);
EXYXY <- Medie2(XY, 2);

VarX <- EXX - EX * EX;
print(VarX);
VarY <- EYY - EY * EY;
print(VarY);
VarXY <- EXYXY - EXY * EXY;
print(VarXY);

Mom1X <- Medie(X, 1);
print(Mom1X);
Mom2X <- Medie(X, 2);
print(Mom2X);
Mom3X <- Medie(X, 3);
print(Mom3X);
Mom4X <- Medie(X, 4);
print(Mom4X);

Mom1Y <- Medie(Y, 1);
print(Mom1Y);
Mom2Y <- Medie(Y, 2);
print(Mom2Y);
Mom3Y <- Medie(Y, 3);
print(Mom3Y);
Mom4Y <- Medie(Y, 4);
print(Mom4Y);

Mom1XY <- Medie2(XY, 1);
print(Mom1XY);
Mom2XY <- Medie2(XY, 2);
print(Mom2XY);
Mom3XY <- Medie2(XY, 3);
print(Mom3XY);
Mom4XY <- Medie2(XY, 4);
print(Mom4XY);

# Computing and checking if exists for the centered moments

MomCentr1X <- MomCentr(X, 1);
if (abs(MomCentr1X) > 1e-8)
{
	print(MomCentr1X);
} else
{
	print("Nu exista");
}
MomCentr2X <- MomCentr(X, 2);
if (abs(MomCentr2X) > 1e-8)
{
	print(MomCentr2X);
} else
{
	print("Nu exista");
}
MomCentr3X <- MomCentr(X, 3);
if (abs(MomCentr3X) > 1e-8)
{
	print(MomCentr3X);
} else
{
	print("Nu exista");
}
MomCentr4X <- MomCentr(X, 4);
if (abs(MomCentr4X) > 1e-8)
{
	print(MomCentr4X);
} else
{
	print("Nu exista");
}

MomCentr1Y <- MomCentr(Y, 1);
if (abs(MomCentr1Y) > 1e-8)
{
	print(MomCentr1Y);
} else
{
	print("Nu exista");
}
MomCentr2Y <- MomCentr(Y, 2);
if (abs(MomCentr2Y) > 1e-8)
{
	print(MomCentr2Y);
} else
{
	print("Nu exista");
}
MomCentr3Y <- MomCentr(Y, 3);
if (abs(MomCentr3Y) > 1e-8)
{
	print(MomCentr3Y);
} else
{
	print("Nu exista");
}
MomCentr4Y <- MomCentr(Y, 4);
if (abs(MomCentr4Y) > 1e-8)
{
	print(MomCentr4Y);
} else
{
	print("Nu exista");
}

MomCentr1XY <- MomCentr2(XY, 1);
if (abs(MomCentr1XY) > 1e-8)
{
	print(MomCentr1XY);
} else
{
	print("Nu exista");
}
MomCentr2XY <- MomCentr2(XY, 2);
if (abs(MomCentr2XY) > 1e-8)
{
	print(MomCentr2XY);
} else
{
	print("Nu exista");
}
MomCentr3XY <- MomCentr2(XY, 3);
if (abs(MomCentr3XY) > 1e-8)
{
	print(MomCentr3XY);
} else
{
	print("Nu exista");
}
MomCentr4XY <- MomCentr2(XY, 4);
if (abs(MomCentr4XY) > 1e-8)
{
	print(MomCentr4XY);
} else
{
	print("Nu exista");
}

# h

# User defined transformation for the random variable

g <- function(x)
{
	return (x + 10);
};

# Computing the E[g(X)] and Var(g(X))

EgX <- Medie(X, 1, g);
print(EgX);
EgXgX <- Medie(X, 2, g);
VargX <- EgXgX - EgX * EgX;
print(VargX);

# User defined transformation for the 2D random variable

g2 <- function(x, y)
{
	return (x + y);
};

# Computing the E[g2(XY)] and Var(g2(XY))

Eg2XY <- Medie2(XY, 1, g2);
print(Eg2XY);
Eg2XYg2XY <- Medie2(XY, 2, g2);
Varg2XY <- Eg2XYg2XY - Eg2XY * Eg2XY;
print(Varg2XY);

# i

# The probability of X in the interval [a; b]

P <- function(X, a, b)
{
	Sum <- 0;

	Step <- 0.1;

	# Integrating

	while (a < b)
	{
		Sum <- Sum + X(a) * Step;

		a <- a + Step;
	}

	return (Sum);
};

# The probability of X in the interval [a; b] and Y in the interval [c; d]

P2 <- function(XY, a, b, c, d)
{
	Step <- 0.1;

	Sum <- 0;

	# Integrating over x

	while (a < b)
	{
		c_copy <- c;

		# Integrating over y

		while (c_copy < d)
		{
			Sum <- Sum + XY(a, c_copy) * Step * Step;

			c_copy <- c_copy + Step;
		}

		a <- a + Step;
	}

	return (Sum);
};

# Testing the P and P2 functions

print(P(X, 0, 1));

print(P2(XY, 0, 1, 0, 1));

# j

# Computing the covariance of X and Y reusing the expected values from g)

CovXY <- EXY - EX * EY;
print(CovXY);

# Computing the corelance coefficient of X and Y reusing the expected values from g)

CoefCorel <- CovXY / sqrt(VarX * VarY);
print(CoefCorel);
