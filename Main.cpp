// Example program
#include<iostream>
#include<string>
#include<cstdlib>
#include<fstream>
double square(double num);
double squareRoot(double num);
double arccos(double num);
double arcsin(double num);
double arctan(double num);
double cos(double num);
double sin(double num);
double tan(double num);
int main()
{
	double pi = 3.14159;
	//double gravity = ;
	//double angle = ;
	//double initialVelocity = ;
        double initialVelocityX = initialVelocity*sin(angle);
        double initialVelocityY = initialVelocity*cos(angle);
	//double timeTillLanded = (-initialVelocity + ) / gravity;
	int numSteps = 9999;
	double* slopes = new double[numSteps];
	double y;
	double x;
	double slope;
	double* ys = new double[numSteps];
	double* xs = new double[numSteps];
        double* vs = new double[numSteps];
        double* angles = new double[numSteps];
	double deltaX = 0.01;
	double angle = pi / 4;
	slope = sin(angle) / cos(angle);
	int count = 0;
	xs[0] = deltaX;
	for (int i = 1;; i++){
		ys[i] = slope *  xs[i]
                if(ys[i]==0){
                    break;
                }
                else if(numSteps+1==10000){
                    numSteps += 99999;
                }
                initialVelocityY -= gravity*(deltaX/initialVelocityY);
                vs[i] = initialVelocityY;
                if(initialVelocity > 0){
                    angle = arctan(initialVelocityY/initialVelocityX);
                }
                else(
                    angle = -arctan(-initialVelocityY/initialVelocityX);
                )
                angles[i] = angle;
		slope = -0.5*(9.81)*square(xs[i] / initialVelocity*cos(angle)) + xs[i] * tan(angle);
		slopes[i] = slope;
		xs[i + 1] -= deltaX;
	}
	fstream file;
	file.open("Output.ext");
	ofstream fout;
	for (int i = 0; i < numSteps + 2; i++){
		fout << xs[i] << " " << ys[i] << endl;
	}
	fstream.close();
	delete[] slopes;
}
double square(double num){
	return num*num;
}
double squareRoot(double num){
	double eps = 0.000001;
	double yo;
	double xo;
	double numYo;
	double numXo;
	double slope;
	if (square(0) != null){
		yo = 0;
		x0 = 0;
	}
	else{
		numYo = yo + eps;
		numXo = xo + eps;
		while (square(num) == null){
			numYo += eps;
		}
		numXo = square(numYo);
	}
	while (abs(numXo) > eps){
		slope = 2 * numXo;
		//(0-numXo) = slope*(numYo - )
		numYo = (0 - numXo) / slope + numYo;
		numXo = square(numYo);
	}
	return numYo;
}
double arccos(double num){
	if (!(-pi / 2 < num && num < pi / 2)){
		cout << "not in the domain of acos" << endl;
		system("pause");
		return null;
	}
	double delta = 0.0000001;
	double integral = 0;
	if (num > 0){
		for (double i = 1 - delta; i >= num + delta; i -= delta){
			integral += ((-1 / (sqrt(1 - pow(i, 2)))) + -1 / (sqrt(1 - pow(i - delta, 2)))) * delta / 2;
		}
		integral = -integral;
	}
	else if (num < 0){
		for (double i = num; i >= -1; i -= delta){
			integral += ((-1 / (sqrt(1 - pow(i, 2)))) + -1 / (sqrt(1 - pow(i - delta, 2)))) * delta / 2;
		}
	}
	else{
		integral = 0;
	}
	return integral;
}
double arcsin(double num){
	if (!(0 < num && num < pi)){
		cout << "not in the domain of acos" << endl;
		system("pause");
		return null;
	}
	double delta = 0.0000001;
	double integral = 0;
	for (double i = 0 + delta; i <= 1; i += delta){
		integral += ((1 / (sqrt(1 - pow(i, 2)))) + 1 / (sqrt(1 - pow(i - delta, 2)))) * delta / 2;
	}
	return integral;
}
double arctan(double num){
	double delta = 0.01;
	double x0;
	double x1;
	double integral = 0;
	double z1;
	double z2;
	bool first = true;
	int count = 0;
	int countCancel = 0;
	int countDone = 0;
	for (int i = -10; i<10; i += 1){
		x0 = i;
		x1 = x0 + delta;
		count = 0;
		first = true;
		integral = 0;
		countDone = 0;
		while (first || abs(integral - 0)>0.000001){
			first = false;
			z1 = sqrt(pow(x0, 2) + pow(y, 2));
			z2 = sqrt(pow(x1, 2) + pow(y, 2));
			integral += -(x0 / sqrt(pow(z1, 2) - pow(x0, 2)) + (x1 / sqrt(pow(z2, 2) - pow(x1, 2)))) * delta / 2;
			if (abs(integral) < 0.01 && count != 0){
				count++;
				if (count>10){
					cout << "Integral " << integral << endl;
					break;
				}
			}
			else{
				countCancel++;
				if (countCancel > 5){
					count = 0;
					countCancel = 0;
					countDone++;
					if (countDone > 2){
						break;
					}
				}
			}
			x0 += delta;
			x1 = x0 + delta;
		}
		cout << integral << " " << x0 << endl;
	}
	return x1;
}
double cos(double num){
	double eps = 0.000001;
	double delta = 0.01;
	if (-pi / 2 < num && num < pi / 2){
		for (double i = -1; i <= 1; i += delta){
			if (abs(arccos(i)) < eps){
				break;
			}
		}
	}
	else{
		cout << "Not in the domain of cos" << endl;
		system("pause");
		return null;
	}
	return i;
        y=cos(x), x=arccos(y)
        f=arccos(y)-x
        df/dy = no more
}
double sin(double num){
	double eps = 0.000001;
	double delta = 0.01;
	if (0 < num && num < pi){
		for (double i = -1; i <= 1; i += delta){
			if (abs(arcsin(i)) < eps){
				break;
			}
		}
	}
	else{
		cout << "Not in the domain of sin" << endl;
		system("pause");
		return null;
	}
	return i;
        
}
double tan(double num){
	double eps = 0.000001;
	double delta = 0.01;
	if (-pi / 2 < num && num < pi / 2){
		for (double i = -1; i <= 1; i += delta){
			if (abs(arctan(i)) < eps){
				break;
			}
		}
	}
	else{
		cout << "Not in the domain of sin" << endl;
		system("pause");
		return null;
	}
	return i;
}
