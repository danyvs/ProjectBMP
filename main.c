#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

struct Pixel {
    unsigned char B, G, R;
};

struct Image {
    unsigned char* header;
    struct Pixel* content;
    unsigned int width, height;
    unsigned int padding;
};

struct Window {
    unsigned int startLine, startColumn;
    unsigned int height, width;
    double correlation;
    struct Pixel color;
};


// ENCRYPTION & DECRYPTION
/**
 *  Generate unsigned integer number with a pseudo random character
 * @param cntNumbers - the number of numbers to be generated
 * @param seed - the number to start the generation with
 * @return a dynamic allocated array with cntNumbers + 1 pseudo random numbers
 */
unsigned int* xorshift32(unsigned int cntNumbers, unsigned int seed) {
    unsigned int* randomNumbers = (unsigned int*)malloc((cntNumbers + 1) * sizeof(unsigned int));
    unsigned int r, i;
    randomNumbers[0] = r = seed;
    for (i = 1; i <= cntNumbers; ++i) {
        r ^= r << 13u;
        r ^= r >> 17u;
        r ^= r << 5u;
        randomNumbers[i] = r;
    }
    return randomNumbers;
}

/**
 *  Load BMP image into internal memory in linearized form
 * @param imagePath - the path of the BMP image
 * @return an Image, a structure who stores information about the image
 */
struct Image loadImageIntoMemory(char *imagePath) {
    struct Image image;
    int i, j;

    FILE* fin = fopen(imagePath, "rb");
    if (fin == NULL) {
        printf("Image was not properly opened!\n");
        exit(EXIT_FAILURE);
    }

    // store information about the image
    image.header = (unsigned char*)malloc(54 * sizeof(unsigned char));
    fread(image.header, sizeof(unsigned char), 54, fin);

    fseek(fin, 18, SEEK_SET);
    fread(&image.width, sizeof(unsigned int), 1, fin);
    fread(&image.height, sizeof(unsigned int), 1, fin);

    image.padding = (image.width % 4) ? 4 - (3 * image.width) % 4 : 0;

    // move back to the beginning of the image
    fseek(fin, 54, SEEK_SET);

    // read the content of the image, starting with the bottom left corner (this is the way BMP images are stored)
    image.content = (struct Pixel*)malloc(image.width * image.height * sizeof(struct Pixel));
    for (i = (int)image.height - 1; i >= 0; --i) {
        // read pixel by pixel, skipping over the padding
        for (j = 0; j < image.width; ++j) {
            struct Pixel pixel;
            fread(&pixel, sizeof(unsigned char), 3, fin);
            image.content[i * image.width + j] = pixel;
        }
        fseek(fin, image.padding, SEEK_CUR);
    }

    fclose(fin);

    return image;
}

/**
 *  Save a BMP image into the external memory (a file)
 * @param filePath - the path where the image is saved
 * @param image - the Image in a linearized form
 */
void saveImageIntoFile(char* filePath, struct Image image) {
    int i;
    unsigned char x = 0;
    FILE* fout = fopen(filePath, "wb");

    fwrite(image.header, sizeof(unsigned char), 54, fout);
    // write the content, starting with the last line and also write the padding
    for (i = (int)image.height - 1; i >= 0; --i) {
        fwrite(&image.content[i * image.width], sizeof(struct Pixel), image.width, fout);
        fwrite(&x, sizeof(unsigned char), image.padding, fout);
    }

    fclose(fout);
}

/**
 *  Using the Durstenfeld algorithm and the pseudo random numbers generated previously, generate a random permutation
 * @param cntNumbers - the numbers of elements of the permutation
 * @param randomNumbers - an array containing the pseudo random numbers generated previously
 * @return an array containing the generated permutation
 */
unsigned int* durstenfeld(unsigned int cntNumbers, const unsigned int* randomNumbers) {
    unsigned int i;

    unsigned int* shuffleArray = (unsigned int*)malloc(cntNumbers * sizeof(unsigned int));
    for (i = 0; i < cntNumbers; ++i)
        shuffleArray[i] = i;

    for (i = cntNumbers - 1; i >= 1; --i) {
        unsigned int randNr = randomNumbers[cntNumbers - i] % (i + 1);
        unsigned int temp = shuffleArray[randNr];
        shuffleArray[randNr] = shuffleArray[i];
        shuffleArray[i] = temp;
    }

    return shuffleArray;
}

/**
 *  Permute the pixels of an image according to given permutation
 * @param image - the Image
 * @param shuffleArray - an array containing the permutation
 */
void shufflePixels(struct Image image, const unsigned int* shuffleArray) {
    unsigned int i;

    // because we may need a permuted pixel, create a temporary array in which to make changes
    struct Pixel* tempImage = (struct Pixel*)malloc(image.height * image.width * sizeof(struct Pixel));
    for (i = 0; i < image.height * image.width; ++i)
        tempImage[shuffleArray[i]] = image.content[i];

    // make the changes in the original image
    for (i = 0; i < image.height * image.width; ++i)
        image.content[i] = tempImage[i];
}

/**
 *  Apply encryption algorithm
 * @param sourceFilePath - filepath of the initial image
 * @param destinationFilePath - filepath where the encrypted image would be stored
 * @param keysFilePath - filepath of the document containing the secret key
 */
void encryptImage(char* sourceFilePath, char* destinationFilePath, char* keysFilePath) {
    FILE* finKeys = fopen(keysFilePath, "r");
    if (!finKeys) {
        printf("File could not be opened!\n");
        return;
    }

    unsigned int randomNumber0, startingValue;
    fscanf(finKeys, "%u%u", &randomNumber0, &startingValue);

    fclose(finKeys);

    struct Image image = loadImageIntoMemory(sourceFilePath);

    unsigned int* randomNumbers = xorshift32(2 * image.height * image.width - 1, randomNumber0);
    unsigned int* shuffleArray = durstenfeld(image.height * image.width, randomNumbers);
    shufflePixels(image, shuffleArray);

    // apply encryption
    image.content[0].R ^= ((startingValue >> (8 * 2)) & 0xFF) ^ ((randomNumbers[image.height * image.width] >> (8 * 2)) & 0xFF);
    image.content[0].G ^= ((startingValue >> (8 * 1)) & 0xFF) ^ ((randomNumbers[image.height * image.width] >> (8 * 1)) & 0xFF);
    image.content[0].B ^= ((startingValue >> (8 * 0)) & 0xFF) ^ ((randomNumbers[image.height * image.width] >> (8 * 0)) & 0xFF);

    unsigned int i;
    for (i = 1; i <= image.height* image.width - 1; ++i) {
        image.content[i].R ^= image.content[i - 1].R ^ ((randomNumbers[image.height * image.width + i] >> (8 * 2)) & 0xFF);
        image.content[i].G ^= image.content[i - 1].G ^ ((randomNumbers[image.height * image.width + i] >> (8 * 1)) & 0xFF);
        image.content[i].B ^= image.content[i - 1].B ^ ((randomNumbers[image.height * image.width + i] >> (8 * 0)) & 0xFF);
    }

    saveImageIntoFile(destinationFilePath, image);

    free(randomNumbers);
    free(shuffleArray);
    free(image.header);
    free(image.content);
}

/**
 *  Calculate the inverse of a permutation
 * @param cntElements - unsigned int, the number of elements of the permutation
 * @param shuffleArray - array of unsigned int, the permutation
 * @return an array of unsigned int, the inverse of the permutation
 */
unsigned int* inversePermutation(unsigned int cntElements, const unsigned int* shuffleArray) {
    unsigned int* tempShuffleArray = (unsigned int*)malloc(cntElements * sizeof(unsigned int));

    unsigned int i;
    for (i = 0; i < cntElements; ++i)
        tempShuffleArray[shuffleArray[i]] = i;

    return tempShuffleArray;
}

/**
 *  Apply decryption algorithm
 * @param filePath - filepath of the initial image
 * @param encryptedFilePath - filepath of the encrypted image
 * @param keysFilePath - filepath of the document containing the secret key
 */
void decryptImage(char* filePath, char* encryptedFilePath, char* keysFilePath) {
    FILE* finKeys = fopen(keysFilePath, "r");
    if (!finKeys) {
        printf("File could not be opened!\n");
        return;
    }

    unsigned int randomNumber0, startingValue;
    fscanf(finKeys, "%u%u", &randomNumber0, &startingValue);

    fclose(finKeys);

    struct Image image = loadImageIntoMemory(encryptedFilePath);

    unsigned int* randomNumbers = xorshift32(2 * image.height * image.width - 1, randomNumber0);
    unsigned int* shuffleArray = durstenfeld(image.height * image.width, randomNumbers);
    shuffleArray = inversePermutation(image.height * image.width, shuffleArray);

    // apply decryption
    unsigned int i;
    for (i = image.height* image.width - 1; i; --i) {
        image.content[i].R ^= image.content[i - 1].R ^ ((randomNumbers[image.height * image.width + i] >> (8 * 2)) & 0xFF);
        image.content[i].G ^= image.content[i - 1].G ^ ((randomNumbers[image.height * image.width + i] >> (8 * 1)) & 0xFF);
        image.content[i].B ^= image.content[i - 1].B ^ ((randomNumbers[image.height * image.width + i] >> (8 * 0)) & 0xFF);
    }

    image.content[0].R ^= ((startingValue >> (8 * 2)) & 0xFF) ^ ((randomNumbers[image.height * image.width] >> (8 * 2)) & 0xFF);
    image.content[0].G ^= ((startingValue >> (8 * 1)) & 0xFF) ^ ((randomNumbers[image.height * image.width] >> (8 * 1)) & 0xFF);
    image.content[0].B ^= ((startingValue >> (8 * 0)) & 0xFF) ^ ((randomNumbers[image.height * image.width] >> (8 * 0)) & 0xFF);

    shufflePixels(image, shuffleArray);

    saveImageIntoFile(filePath, image);

    free(randomNumbers);
    free(shuffleArray);
    free(image.header);
    free(image.content);
}


/**
 *  Print the values of the chi-squared test for a BMP image on every color channel (R, G, B)
 * @param imagePath - the path of the image
 */
void printChiSquareTest(char* imagePath) {
    struct Image image = loadImageIntoMemory(imagePath);

    // theoretically estimated frequency
    double fbar = image.height * image.width / 256.0;

    // chi-squared test value for every channel
    double chiR = 0, chiG = 0, chiB = 0;

    unsigned int i, j, intensity;
    for (intensity = 0; intensity < 256; ++intensity) {
        double sumR = 0, sumG = 0, sumB = 0;
        for (i = 0; i < image.height; ++i)
            for (j = 0; j < image.width; ++j) {
                // update frequencies
                if (image.content[i * image.width + j].R == intensity)
                    ++sumR;

                if (image.content[i * image.width + j].G == intensity)
                    ++ sumG;

                if (image.content[i * image.width + j].B == intensity)
                    ++sumB;
            }

        chiR += (sumR - fbar) * (sumR - fbar) / fbar;
        chiG += (sumG - fbar) * (sumG - fbar) / fbar;
        chiB += (sumB - fbar) * (sumB - fbar) / fbar;
    }

    printf("R: %0.2lf\n", chiR);
    printf("G: %0.2lf\n", chiG);
    printf("B: %0.2lf\n", chiB);

    free(image.header);
    free(image.content);
}


// TEMPLATE MATCHING
/**
 *  Transform color image into grayscale image
 * @param sourceFilePath - filepath of the source image
 * @param destinationFilePath - filepath of the grayscale image
 */
void grayscaleImage(char* sourceFilePath, char* destinationFilePath) {
    struct Image image = loadImageIntoMemory(sourceFilePath);

    unsigned int i, j;
    for (i = 0; i < image.height; ++i)
        for (j = 0; j < image.width; ++j) {
            unsigned int idx = i * image.width + j;
            unsigned char gray = 0.299 * image.content[idx].R + 0.587 * image.content[idx].G + 0.114 * image.content[i].B;
            image.content[idx].R = image.content[idx].G = image.content[idx].B = gray;
        }

    saveImageIntoFile(destinationFilePath, image);

    free(image.header);
    free(image.content);
}

/**
 *  Calculate correlation between 2 images with a given formula
 * @param image - Image, representing the original mage
 * @param template - Image, representing the template
 * @param line - int, the beginning line of the window
 * @param column - int, the beginning column of the window
 * @return correlation, a double
 */
double calculateCorrelation(struct Image image, struct Image template, unsigned int line, unsigned int column) {
    unsigned int i, j;
    unsigned int cntPixels = template.height * template.width;

    // mean value of pixel's grayscale intensity in the window
    double meanPixelsIntensityImage = 0;
    for (i = 0; i < template.height; ++i)
        for (j = 0; j < template.width; ++j)
            meanPixelsIntensityImage += image.content[(i + line) * image.width + (j + column)].R;
    meanPixelsIntensityImage /= cntPixels;

    // mean value of pixel's grayscale intensity in the template
    double meanPixelsIntensityTemplate = 0;
    for (i = 0; i < template.height; ++i)
        for (j = 0; j < template.width; ++j)
            meanPixelsIntensityTemplate += template.content[i * template.width + j].R;
    meanPixelsIntensityTemplate /= cntPixels;

    // standard deviation of pixel's grayscale intensity in the window
    double standardDeviationImage = 0;
    for (i = 0; i < template.height; ++i)
        for (j = 0; j < template.width; ++j) {
            double pixelDeviation = image.content[(i + line) * image.width + (j + column)].R - meanPixelsIntensityImage;
            pixelDeviation *= pixelDeviation;
            standardDeviationImage += pixelDeviation;
        }
    standardDeviationImage /= (cntPixels - 1);
    standardDeviationImage = sqrt(standardDeviationImage);

    // standard deviation of pixel's grayscale intensity in the template
    double standardDeviationTemplate = 0;
    for (i = 0; i < template.height; ++i)
        for (j = 0; j < template.width; ++j) {
            double pixelDeviation = template.content[i * template.width + j].R - meanPixelsIntensityTemplate;
            pixelDeviation *= pixelDeviation;
            standardDeviationTemplate += pixelDeviation;
        }
    standardDeviationTemplate /= (cntPixels - 1);
    standardDeviationTemplate = sqrt(standardDeviationTemplate);

    double correlation = 0;
    for (i = 0; i < template.height; ++i)
        for (j = 0; j < template.width; ++j) {
            unsigned char intensityPixelImage = image.content[(i + line) * image.width + (j + column)].R;
            unsigned char intensityPixelTemplate = template.content[i * template.width + j].R;

            double currentCorrelation = intensityPixelImage - meanPixelsIntensityImage;
            currentCorrelation *= (intensityPixelTemplate - meanPixelsIntensityTemplate);
            currentCorrelation /= standardDeviationImage;
            currentCorrelation /= standardDeviationTemplate;
            correlation += currentCorrelation;
        }
    correlation /= cntPixels;

    return correlation;
}

/**
 *  Get all the matchings for a given template who have a correlation higher than the threshold
 * @param imagePath - the filepath of the initial image
 * @param templatePath - the filepath of the template
 * @param threshold - double, the threshold for a matching to be considered correct
 * @param matches - array of Window, all the found matchings for the given image and the given template
 * @param cntMatches - unsigned int, the number of matchings
 */
void templateMatching(char* imagePath, char* templatePath, double threshold, struct Window** matches, unsigned int* cntMatches) {
    struct Image image = loadImageIntoMemory(imagePath);
    struct Image template = loadImageIntoMemory(templatePath);

    *matches = NULL;
    *cntMatches = 0;

    unsigned int i, j;
    for (i = 0; i + template.height < image.height; ++i)
        for (j = 0; j + template.width < image.width; ++j) {
            double corr = calculateCorrelation(image, template, i, j);
            if (corr >= threshold) {
                ++(*cntMatches);
                struct Window* tempWindows = (struct Window*)realloc((*matches), (*cntMatches) * sizeof(struct Window));
                if (!tempWindows) {
                    printf("Could not relocate memory!\n");
                    free(*matches);
                    exit(EXIT_FAILURE);
                }
                else {
                    *matches = tempWindows;
                    (*matches)[*cntMatches - 1].startLine = i;
                    (*matches)[*cntMatches - 1].startColumn = j;
                    (*matches)[*cntMatches - 1].height = template.height;
                    (*matches)[*cntMatches - 1].width = template.width;
                    (*matches)[*cntMatches - 1].correlation = corr;
                }
            }
        }

    free(image.header);
    free(image.content);
    free(template.header);
    free(template.content);
}

/**
 *  Draw border of a window with a given color
 * @param image - the Image
 * @param window - the Window to be bordered
 * @param color - a Pixel, representing the color
 */
void drawBorderWindow(struct Image image, struct Window window, struct Pixel color) {
    unsigned int i;

    // vertically
    for (i = 0; i < window.height; ++i) {
        image.content[(i + window.startLine) * image.width + window.startColumn] = color;
        image.content[(i + window.startLine) * image.width + (window.startColumn + window.width - 1)] = color;
    }
    // horizontally
    for (i = 0; i < window.width; ++i) {
        image.content[window.startLine * image.width + window.startColumn + i] = color;
        image.content[(window.startLine + window.height - 1) * image.width + window.startColumn + i] = color;
    }
}

/**
 *  Generate an array of Pixels with the color to be bordered each number
 * @return array of Pixels
 */
struct Pixel* initColorsForPixels() {
    struct Pixel* colorOfNumbers = (struct Pixel*)malloc(10 * sizeof(struct Pixel));

    colorOfNumbers[0].R = 255, colorOfNumbers[0].G =   0, colorOfNumbers[0].B =   0;
    colorOfNumbers[1].R = 255, colorOfNumbers[1].G = 255, colorOfNumbers[1].B =   0;
    colorOfNumbers[2].R =   0, colorOfNumbers[2].G = 255, colorOfNumbers[2].B =   0;
    colorOfNumbers[3].R =   0, colorOfNumbers[3].G = 255, colorOfNumbers[3].B = 255;
    colorOfNumbers[4].R = 255, colorOfNumbers[4].G =   0, colorOfNumbers[4].B = 255;
    colorOfNumbers[5].R =   0, colorOfNumbers[5].G =   0, colorOfNumbers[5].B = 255;
    colorOfNumbers[6].R = 192, colorOfNumbers[6].G = 192, colorOfNumbers[6].B = 192;
    colorOfNumbers[7].R = 255, colorOfNumbers[7].G = 140, colorOfNumbers[7].B =   0;
    colorOfNumbers[8].R = 128, colorOfNumbers[8].G =   0, colorOfNumbers[8].B = 128;
    colorOfNumbers[9].R = 128, colorOfNumbers[9].G =   0, colorOfNumbers[9].B =   0;

    return colorOfNumbers;
}

/**
 *  Get all the matches that have a matching of at least an input threshold
 * @param allMatches - array of Window, stores the matches
 * @param cntAllMatches - int, stores the number of matches
 */
void getAllMatches(struct Window** allMatches, unsigned int *cntAllMatches) {
    struct Pixel* colorOfNumbers = initColorsForPixels();

    char filePath[20];
    printf("Input the filepath of the initial image: ");
    scanf("%s", filePath);

    FILE* fin = fopen(filePath, "r");
    if (!fin) {
        printf("File could not be opened!\n");
        exit(EXIT_FAILURE);
    }

    // make the image grayscale
    fscanf(fin, "%s", filePath);
    grayscaleImage(filePath, "grayscaleImage.bmp");

    // input the threshold
    double threshold;
    printf("Input the threshold: ");
    scanf("%lf", &threshold);

    // read the number of templates
    unsigned int cntTemplates;
    fscanf(fin, "%u", &cntTemplates);

    *allMatches = NULL;
    *cntAllMatches = 0;

    unsigned int i, j;
    for (i = 0; i < cntTemplates; ++i) {
        // read the template and make it grayscale
        fscanf(fin, "%s", filePath);
        grayscaleImage(filePath, "grayscaleTemplate.bmp");

        // find all the matches for the current template and store them in a temporary array
        unsigned int cntMatches;
        struct Window* matches;
        templateMatching("grayscaleImage.bmp", "grayscaleTemplate.bmp", threshold, &matches, &cntMatches);

        // store the color for current matches
        for (j = 0; j < cntMatches; ++j)
            matches[j].color = colorOfNumbers[i];

        // move the matchings for the current template from the temporary array to the permanent array
        *cntAllMatches += cntMatches;
        struct Window* tempAllMatches = (struct Window*)realloc(*allMatches, (*cntAllMatches) * sizeof(struct Window));
        if (!tempAllMatches && cntMatches) {
            printf("Could not relocate memory!\n");
            free(*allMatches);
            exit(EXIT_FAILURE);
        }
        else {
            *allMatches = tempAllMatches;
            for (j = 0; j < cntMatches; ++j)
                (*allMatches)[*cntAllMatches - cntMatches + j] = matches[j];
        }
    }

    fclose(fin);
}

unsigned int min(unsigned int a, unsigned int b) {
    return a < b ? a : b;
}

unsigned int max(unsigned int a, unsigned int b) {
    return a > b ? a : b;
}

/**
 *  Helper function to calculate the intersection of 2 windows
 */
double intersection(struct Window a, struct Window b) {
    // (x1, y1) - left top corner
    unsigned int x1 = max(a.startLine, b.startLine);
    unsigned int y1 = max(a.startColumn, b.startColumn);
    // (x2, y2) - right bottom corner
    unsigned int x2 = min(a.startLine + a.height, b.startLine + b.height);
    unsigned int y2 = min(a.startColumn + a.width, b.startColumn + b.width);

    double res = 0;
    if (a.startLine <= x1 && x1 <= a.startLine + a.height && a.startColumn <= y1 && y1 <= a.startColumn + a.width)
        if (a.startLine <= x2 && x2 <= a.startLine + a.height && a.startColumn <= y2 && y2 <= a.startColumn + a.width) {
            res = (x2 - x1) * (y2 - y1);
            res /= (a.height * a.width + b.height * b.width - res);
        }

    return res;
}

/**
 *  Helper function for qsort, sorting by correlation
 */
int compareByCorrelation(const void* a, const void* b) {
    struct Window x = *(struct Window*)a;
    struct Window y = *(struct Window*)b;
    if (x.correlation < y.correlation)
        return 1;
    else
        return -1;
}

/**
 *  Keep the matches with the highest correlation and remove the others
 * @param allMatches - array of Window, stores all the matches
 * @param cntAllMatches - unsigned int, the number of matches
 */
void nonMaximalElimination(struct Window** allMatches, unsigned int* cntAllMatches) {
    unsigned int i, j;

    // to remove elements efficiently, create an array who stores if a match should be removed (1) or not (0)
    int* eliminated = (int*)malloc(*cntAllMatches * sizeof(int));
    for (i = 0; i < *cntAllMatches; ++i)
        eliminated[i] = 0;

    for (i = 0; i < *cntAllMatches; ++i)
        if (!eliminated[i])
            for (j = i + 1; j < *cntAllMatches; ++j)
                if (!eliminated[j])
                    if (intersection((*allMatches)[i], (*allMatches)[j]) > 0.2)
                        eliminated[j] = 1;

    // keep only the "best" matches
    unsigned int index = 0;
    for (i = 0; i < (*cntAllMatches); ++i)
        if (!eliminated[i])
            (*allMatches)[index++] = (*allMatches)[i];
    *cntAllMatches = index;

    struct Window* tempAllMatches = (struct Window*)realloc(*allMatches, (*cntAllMatches) * sizeof(struct Window));
    if (!tempAllMatches) {
        printf("Could not relocate memory!\n");
        free(*allMatches);
        exit(EXIT_FAILURE);
    }
    else
        *allMatches = tempAllMatches;
}

/**
 *  Draw borders of the windows
 * @param sourcefilePath - initial image's filepath
 * @param destinationFilePath - saved image's filepath
 * @param allMatches - array, containing all the matches
 * @param cntAllMatches - unsigned int, the number of windows matching
 */
void drawBorders(char* sourcefilePath, char* destinationFilePath, struct Window* allMatches, unsigned int cntAllMatches) {
    struct Image image = loadImageIntoMemory(sourcefilePath);

    unsigned int i;
    for (i = 0; i < cntAllMatches; ++i)
        drawBorderWindow(image, allMatches[i], allMatches[i].color);

    saveImageIntoFile(destinationFilePath, image);
}

int main() {
    char imagePath[20], encryptedImagePath[20], secretKeysTextPath[20], finalImagePath[20];

    // ENCRYPTION & DECRYPTION
    printf("Input filepath of the initial image: ");
    scanf("%s", imagePath);
    printf("Input filepath of the encrypted image: ");
    scanf("%s", encryptedImagePath);
    printf("Input filepath of the secret key: ");
    scanf("%s", secretKeysTextPath);

    encryptImage(imagePath, encryptedImagePath, secretKeysTextPath);
    decryptImage(imagePath, encryptedImagePath, secretKeysTextPath);

    printf("The values of the chi-squared test for the initial image are:\n");
    printChiSquareTest(imagePath);
    printf("The values of the chi-squared test for the encrypted image are:\n");
    printChiSquareTest(encryptedImagePath);

    // TEMPLATE MATCHING
    unsigned int cntAllMatches;
    struct Window* allMatches;

    getAllMatches(&allMatches, &cntAllMatches);

    qsort(allMatches, cntAllMatches, sizeof(struct Window), compareByCorrelation);
    nonMaximalElimination(&allMatches, &cntAllMatches);

    printf("Input filepath of the image on which to draw: ");
    scanf("%s", finalImagePath);
    drawBorders(imagePath, finalImagePath, allMatches, cntAllMatches);

    return 0;
}