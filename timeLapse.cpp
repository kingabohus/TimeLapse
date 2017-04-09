// timeLapse.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#ifdef _WIN32
#define _CRT_SECURE_NO_DEPRECATE
#endif
#include <stdio.h>
#pragma warning (disable : 4996)
#include <iostream>
#include <fstream>
#include <opencv2/opencv.hpp>
#include <vector>
#include <math.h>
#include <string>
#include <sstream>

#define PI 3.14159265
#define MATCHTRESHOLD 0.4
#define MAXTILT 4
#define CANNYT 100
#define SMOOTHING 0.8
#define CORRTH 0.5
#define SHRINK 6

using namespace cv;
using namespace std;

double correlation(cv::Mat &image_1, cv::Mat &image_2) {

	// convert data-type to "float"
	cv::Mat im_float_1;
	image_1.convertTo(im_float_1, CV_32F);
	cv::Mat im_float_2;
	image_2.convertTo(im_float_2, CV_32F);

	int n_pixels = im_float_1.rows * im_float_1.cols;

	// Compute mean and standard deviation of both images
	cv::Scalar im1_Mean, im1_Std, im2_Mean, im2_Std;
	meanStdDev(im_float_1, im1_Mean, im1_Std);
	meanStdDev(im_float_2, im2_Mean, im2_Std);

	// Compute covariance and correlation coefficient
	double covar = (im_float_1 - im1_Mean).dot(im_float_2 - im2_Mean) / n_pixels;
	double correl = covar / (im1_Std[0] * im2_Std[0]);

	return correl;
}

void getEdgePoints(const Mat& image, Mat& draw)
{
	//Mat blurred;
	//medianBlur(image, blurred, 5);
	Mat gray, edge;
	cvtColor(image, gray, CV_BGR2GRAY);
	Canny(gray, edge, CANNYT, CANNYT*2, 3);
	edge.convertTo(draw, CV_8U);
}

void rotate(const Mat& src, Mat& dst, double  angle, bool black)
{
	Point2d pt(src.cols / 2., src.rows / 2.);
	Mat r = getRotationMatrix2D(pt, angle, 1.0);
	if (black) { warpAffine(src, dst, r, src.size(), INTER_LINEAR, BORDER_CONSTANT); }
	else { warpAffine(src, dst, r, src.size(), INTER_LINEAR, BORDER_CONSTANT, Scalar(192, 192, 192)); }
}

void translate(const Mat& img, Mat& trans_img, double offsetx, double offsety, bool black)
{
	Mat trans_mat = (Mat_<double>(2, 3) << 1, 0, offsetx, 0, 1, offsety);
	if (black) { warpAffine(img, trans_img, trans_mat, img.size(), INTER_LINEAR, BORDER_CONSTANT); }
	else { warpAffine(img, trans_img, trans_mat, img.size(), INTER_LINEAR, BORDER_CONSTANT, Scalar(192, 192, 192)); }
}

void padToSize(const Mat& img, Mat& trans_img, Size size, bool black)
{
	int h = size.height;
	int w = size.width;
	//if (h > img.rows)
	//{
		int left = (w - img.cols) / 2;
		int right = ((w - img.cols) / 2) + ((w - img.cols) % 2);
		int top = (h - img.rows) / 2;
		int bottom = ((h - img.rows) / 2) + ((h - img.rows) % 2);
		if (black) {copyMakeBorder(img, trans_img, top, bottom, left, right, BORDER_CONSTANT);}
		else { copyMakeBorder(img, trans_img, top, bottom, left, right, BORDER_CONSTANT, Scalar(192, 192, 192)); }
	//}
	//else
	//{
	//	Mat newe2;

	//}
}

class Transformer
{
	Size m_size;
	double m_x;
	double m_y;
	float m_d;
	

public:
	Transformer(Size size, double x = 0, double y = 0, double d = 0)
	{
		m_size = size;
		m_x = x;
		m_y = y;
		m_d = d;
	}

	Mat applyTransform(const Mat& img, double tx, double ty, double angle)
	{
		m_x = (m_x + tx)*0.95;
		m_y = (m_y + ty)*0.95;
		m_d = (m_d + angle)*0.95;

		Mat trans_img;
		Mat rotated;
		Mat padded;
		
		padToSize(img, padded, m_size, false);
		rotate(padded, rotated, m_d, false);
		translate(rotated, trans_img, m_x, m_y, false);
		return trans_img;
	}

	Mat trialtransform(const Mat& img, double tx, double ty)
	{
		Mat trans_img;
		Mat padded;
		padToSize(img, padded, img.size(), false);
		translate(padded, trans_img, m_x + tx, m_y + ty, false);

		return trans_img;
	}

	double getX()
	{
		return m_x;
	}

	double getY()
	{
		return m_y;
	}

	double getAngle()
	{
		return m_d;
	}

	void reset()
	{
		m_x = 0;
		m_y = 0;
		m_d = 0;
	}
};



std::vector<std::string> readFilenames(std::string filename)
{
	ifstream inputFile;

	//Declare other variables (forgot to add these in my previous EDIT, sorry)
	int number_of_files;
	string line;

	//Open file list and count number of files
	inputFile.clear();
	inputFile.open(filename, ios::in);

	//exit and prompt error message if file could not be opened
	if (!inputFile) {
		cerr << "File list could not be opened" << endl;
		exit(1);
	}// end if
	number_of_files = 0;
	while (getline(inputFile, line))
		number_of_files++;
	cout << "Number of pictures: " << number_of_files << endl;
	std::vector<std::string> filelist(number_of_files);
	inputFile.close();
	//Re-open file list and store filenames in a string array
	inputFile.clear();
	inputFile.open(filename, ios::in);
	//exit and prompt error message if file could not be opened
	if (!inputFile) {
		cerr << "File list could not be opened" << endl;
		exit(1);
	}// end if
	 // store filenames
	int i;
	i = 0;
	while (getline(inputFile, line)) {
		filelist[i] = line;
		//cout << filelist[i] << endl;
		i = i + 1;
	}
	inputFile.close();
	return filelist;
}

int main(int argc, char** argv)
{
	if (argc != 2)
	{
		cout << "Argument: file with filenames" << endl;
		return -1;
	}
	//Declarations for I/O files
	std::vector<std::string> filelist = readFilenames(argv[1]);
	Mat img_1, img_2;
	stringstream fileIn;
	fileIn << "./input/" << filelist[0].c_str();
	img_1 = imread(fileIn.str(), CV_LOAD_IMAGE_COLOR);// Read the first pic
	
	Size size = img_1.size();
	Size sizeOut(size.width+100, size.height+100);
	Transformer transformer(sizeOut, 0, 0, 0);
	
	stringstream file;
	file << "./output/trimg0.ppm";
	Mat ti;
	padToSize(img_1, ti, sizeOut, false);
	imwrite(file.str(), ti);
	int w = size.width;
	int h = size.height;
	auto writer = cv::VideoWriter("out.avi", VideoWriter::fourcc('M', 'J', 'P', 'G'), 12, sizeOut);
	writer << ti;
	
	//namedWindow("Image1", WINDOW_AUTOSIZE);
	//imshow("Image1", img_1);

	for (int i = 1; i < filelist.size(); i++)
	//for (int i = 420; i < 430; i++)
	{
		std::cout << "Image #" << i << "\n";
		Mat halfImg_1;
		resize(img_1, halfImg_1, Size(), 1.0/SHRINK, 1.0/SHRINK, INTER_LINEAR);
		Mat edges1;
		getEdgePoints(halfImg_1, edges1);
		int edgeTotal1 = countNonZero(edges1);
		
		stringstream fileIn;
		fileIn << "./input/" << filelist[i].c_str();
		img_2 = imread(fileIn.str(), CV_LOAD_IMAGE_COLOR);
		if (w < img_2.cols)
		{
			Size newSize(w, h*(img_2.cols/w));
			resize(img_2, img_2, newSize);
			
		}
		Mat halfImg_2;
		resize(img_2, halfImg_2, Size(), 1.0 / SHRINK, 1.0 / SHRINK, INTER_LINEAR);
		Size halfSize = halfImg_2.size();
		
		Mat edges2;
		getEdgePoints(halfImg_2, edges2);

		/*if (edges2.size() != edges1.size())
		{
			if (edges1.rows > edges2.rows)
			{

				padToSize(edges2, edges2, edges1.size(), true);
			}
			else
			{
				padToSize(edges1, edges1, edges2.size(), true);
			}
		}*/

		double bestMatch = 0;
		double bestangle = -MAXTILT;
		double bestx = 0;
		double besty = 0;
		Mat hann;
		createHanningWindow(hann, halfSize, CV_32F);

		for (double d = -MAXTILT; d <= MAXTILT; d += 0.1)
		{ 
			Mat rotimg2;
			rotate(edges2, rotimg2, d, true);
			Mat e32f1, e32f2;
			rotimg2.convertTo(e32f2, CV_32FC1);
			edges1.convertTo(e32f1, CV_32FC1);
			
			Point2d pt = phaseCorrelate(e32f1, e32f2, hann);
			Mat rottrans2;
			translate(rotimg2, rottrans2, -pt.x, -pt.y, true);
			Mat both = edges1 & rottrans2;
			int edgeTotal2 = countNonZero(rottrans2);
			int matchNonZero = countNonZero(both);

			double match = max(matchNonZero*1.0/edgeTotal1, matchNonZero*1.0 / edgeTotal1
			)  ;
			if (match > bestMatch)
			{
				bestMatch = match;
				bestangle = d;
				bestx = -pt.x;
				besty = -pt.y;
				
			}
		}

			/*double cx;
			double cy;
			double bestCorr = 0;
			for (int px = -1; px <= 1; px += 1)
			{
				for (int py = -1; py <= 1; py += 1)
				{
					Mat tr2 = transformer.trialTransform(img_2, SHRINK*bestx + px, SHRINK*besty + py);
					double corr = correlation(img_1, tr2);
					std::cout << "Best corr: " << corr << "\n";
					if (corr > bestCorr)
					{
						bestCorr = corr;
						cx = px;
						cy = py;
					}

				}

			}*/
		//Mat trimg2 = transformer.applyTransform(img_2, SHRINK*bestx, SHRINK *besty, bestangle);
		Mat tr2 = transformer.trialtransform(img_2, SHRINK * bestx, SHRINK * besty);
		double corr = correlation(img_1, tr2);
			std::cout << "corr: " << corr << "\n";

			//std::cout << "Best match: " << bestx << ", " << besty << " angle = " << bestangle << " Match = " << bestMatch <<", " << edgeTotal1 << "\n";
			//cout << bestMatch*edgeTotal1 / edgeTotal2 << "\n";

			if (corr > CORRTH)
			{
				std::cout << "Best match: " << SHRINK *bestx << ", " << SHRINK *besty << " angle = " << bestangle <<  "\n";
				Mat trimg2 = transformer.applyTransform(img_2, SHRINK*bestx, SHRINK *besty, bestangle);
				stringstream file;
				file << "./output/trimg" << i << ".ppm";
				imwrite(file.str(), trimg2);
				/*if (trimg2.size() != size)
				{
					trimg2 = padToSize(trimg2, size);
				}*/
				writer << trimg2;
				
			}
			else
			{
				transformer.reset();
				cout << "RESET TO CENTER\n";
				Mat black(sizeOut.height, sizeOut.width, CV_8UC3, Scalar(192, 192, 192));
				writer << black;
				Mat img2out;
				padToSize(img_2, img2out, sizeOut, false);
				stringstream file;
				file << "./output/trimg" << i << ".ppm";
				imwrite(file.str(), img2out);
				//if (img_2.size() != size)
				//{
				//	img_2 = padToSize(img_2, size);
				//}
				writer << img2out;
				
			}
		
		img_1 = img_2;
	}
	writer.release();
	waitKey(0); // Wait for a keystroke in the window
	return 0;
}

