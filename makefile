# Nom de l'exécutable à générer
TARGET = program

# Compilateur
CC = gcc

# Options de compilation
CFLAGS = -Wall 

# Fichiers source
SOURCES = main.c algebre.c fonction.c BiCGstab.c parametre.c

# Fichiers objet générés à partir des fichiers source
OBJECTS = $(SOURCES:.c=.o)

# Règle principale
all: $(TARGET)

# Comment générer l'exécutable
$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(TARGET) -lm

# Règle générique pour transformer un .c en .o
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Nettoyage des fichiers générés
clean:
	rm -f $(OBJECTS) $(TARGET) *.dat
