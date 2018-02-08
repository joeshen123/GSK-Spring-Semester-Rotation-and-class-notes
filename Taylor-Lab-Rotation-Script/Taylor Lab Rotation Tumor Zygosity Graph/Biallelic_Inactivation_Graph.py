import numpy as np
import cv2

# Draw the Germline Biallelic events
image = np.zeros([800,1000,3], dtype = np.uint8)

image.fill(255)
#Draw the Germline chromosome
cv2.line(image, (56,439), (56, 280), (0,0,255), 15)
cv2.line(image, (120,439), (120, 280), (0,0,255), 15)
cv2.line(image, (30,346), (74, 386), (255,0,0), 8)
cv2.line(image, (30,386), (74, 346), (255,0,0), 8)

#Draw the Somatic chromosome one
cv2.line(image, (696,110), (696, 269), (0,0,255), 15)
cv2.line(image, (760,110), (760, 269), (0,0,255), 15)
cv2.line(image, (673,142), (715, 169), (255,0,0), 8)
cv2.line(image, (673,169), (715, 142), (255,0,0), 8)
cv2.line(image, (730,138), (785, 169), (255,0,0), 8)
cv2.line(image, (730,169), (785, 138), (255,0,0), 8)

#Draw the Somatic chromosome two
cv2.line(image, (696,496), (696, 655), (0,0,255), 15)
cv2.line(image, (760,496), (760, 530), (0,0,255), 15)
cv2.line(image, (760,600), (760, 655), (0,0,255), 15)
cv2.line(image, (673,526), (715, 556), (255,0,0), 8)
cv2.line(image, (673,556), (715, 526), (255,0,0), 8)

#Put the arrow line
cv2.arrowedLine(image,(140,363), (560,363), (255,0,0), 6)

#Put the text
font = cv2.FONT_HERSHEY_SIMPLEX
cv2.putText(image, "Germline Mutation", (30,211), font, 1.1, (0,0,0), 4, cv2.LINE_4 )
cv2.putText(image, "Normal Cell", (30,761), font, 1.1, (20,0,100), 4, cv2.LINE_4 )
cv2.putText(image, "Tumor Cell", (620,761), font, 1.1, (20,0,100), 4, cv2.LINE_4 )
cv2.putText(image, "Second Somatic Mutation", (500,66), font, 1.1, (0,0,0), 4, cv2.LINE_4 )
cv2.putText(image, "Loss of heterozygosity", (500,450), font, 1.1, (0,0,0), 4, cv2.LINE_4 )


cv2.imwrite('Germline_Image.tiff', image)


# Draw the Somatic Biallelic events
image_S = np.zeros([800,1200,3], dtype = np.uint8)

image_S.fill(255)

#Draw the Germline chromosome
cv2.line(image_S, (56,439), (56, 280), (0,0,255), 15)
cv2.line(image_S, (120,439), (120, 280), (0,0,255), 15)
cv2.line(image_S, (30,346), (74, 386), (255,0,0), 8)
cv2.line(image_S, (30,386), (74, 346), (255,0,0), 8)

#Draw the Somatic chromosome one
cv2.line(image_S, (696,110), (696, 269), (0,0,255), 15)
cv2.line(image_S, (760,110), (760, 269), (0,0,255), 15)
cv2.line(image_S, (673,142), (715, 169), (255,0,0), 8)
cv2.line(image_S, (673,169), (715, 142), (255,0,0), 8)
cv2.line(image_S, (730,138), (785, 169), (255,0,0), 8)
cv2.line(image_S, (730,169), (785, 138), (255,0,0), 8)

#Draw the Somatic chromosome two
cv2.line(image_S, (696,496), (696, 655), (0,0,255), 15)
cv2.line(image_S, (760,496), (760, 530), (0,0,255), 15)
cv2.line(image_S, (760,600), (760, 655), (0,0,255), 15)
cv2.line(image_S, (673,526), (715, 556), (255,0,0), 8)
cv2.line(image_S, (673,556), (715, 526), (255,0,0), 8)

#Draw the Somatic chromosome two
cv2.line(image_S, (973,308), (973, 350), (0,0,255), 15)
cv2.line(image_S, (973,430), (973, 467), (0,0,255), 15)

cv2.line(image_S, (1037,308), (1037, 350), (0,0,255), 15)
cv2.line(image_S, (1037,430), (1037, 467), (0,0,255), 15)


#Draw the arrow
cv2.arrowedLine(image_S,(140,363), (560,363), (255,0,0), 6)

#Add text
font = cv2.FONT_HERSHEY_SIMPLEX
cv2.putText(image_S, "Somatic Mutation", (30,211), font, 1.1, (0,0,0), 4, cv2.LINE_4 )
cv2.putText(image_S, "Tumor Cell", (30,761), font, 1.1, (20,0,100), 4, cv2.LINE_4 )
cv2.putText(image_S, "Tumor Cell", (620,761), font, 1.1, (20,0,100), 4, cv2.LINE_4 )
cv2.putText(image_S, "Tumor Cell", (950,761), font, 1.1, (20,0,100), 4, cv2.LINE_4 )
cv2.putText(image_S, "Second Mutation hit ", (500,66), font, 1.1, (0,0,0), 4, cv2.LINE_4 )
cv2.putText(image_S, "Loss of heterozygosity", (500,450), font, 1.1, (0,0,0), 4, cv2.LINE_4 )

#use for loop to add newline in text
text = "Homozygous \n Deletion"

y0 = 196

yd = 50

for i,line in enumerate(text.split('\n')):
    cv2.putText(image_S, line, (930, y0 + i * yd), font, 1.1, (0, 0, 0), 4, cv2.LINE_4)

#Add line to seperate Homodel from Others
cv2.line(image_S, (899,0), (899,1200), (0,100,200),5,8)

cv2.imwrite('Somatic_Image.tiff', image_S)

