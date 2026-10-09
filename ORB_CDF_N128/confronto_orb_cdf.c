#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdint.h> //per utilizzare uint64_t, un intero senza segno formato esattamente da 64 bit
//si userà al posto di Erank[N=128], perché, così ragionando, dentro ogni nodo della min-heap sarebbero contenuti, per N=128 => 128*4=512 byte, sono troppi soprattutto per un budget massimo di 2500000 query
//p.e con 500000 candidati si arriva a 260 MB

#define N 128 //lunghezza della codeword
#define NH 64 //numero di righe della matrice di controllo di parità

//RAPPRESENTAZIONE SPARSA DELLA MATRICE H
//per ogni riga memorizza gli indici delle 8 colonne contenenti un 1
int H_pos[NH][8];

//CONTATORE DEL NUMERO DI QUERY EFFETTUATE DA CDF-ORBGRAND
//NUOVO RISPETTO AL VECCHIO PROGETTO:
//ogni query corrisponde a un pattern candidato effettivamente sottoposto al controllo di sindrome
unsigned long long numero_query = 0;

//NUMERO MASSIMO DI QUERY CONSENTITE A CDF-ORBGRAND
//se il decoder raggiunge questo numero di query senza aver trovato una codeword valida, la ricerca viene interrotta e il blocco viene classificato come abbandono

#define MAX_QUERY 2500000ULL

//FLAG PER ABILITARE O DISABILITARE LE STAMPE DETTAGLIATE
//1 = stampa completa utilizzata per il debug e la validazione
//0 = output ridotto per le simulazioni su molti blocchi
#define STAMPE_DETTAGLIATE 0

//RAPPRESENTAZIONE COMPATTA DEI PATTERN CDF PER N=128
//Nel progetto toy ogni pattern Erank veniva memorizzato tramite un vettore di N interi
//Con N=128 questa soluzione richiederebbe molta memoria, soprattutto perché la min-heap può contenere molti nodi
//Un uint64_t contiene esattamente 64 bit => due uint64_t permettono di rappresentare tutti i 128 ranghi usando 16 byte

typedef struct {

	uint64_t Erank_basso; //ranghi di affidabilità da 1 a 64
	uint64_t Erank_alto;  //ranghi di affidabilità da 65 a 128

	double peso_cdf;      //peso complessivo del pattern secondo la metrica CDF
	int rango_massimo;    //rango di affidabilità più grande presente nel pattern

} NodoHeapCDF;


//STRUTTURA CHE RAPPRESENTA LA MIN-HEAP CDF
//nodi è un puntatore all'array di nodi che verrà allocato dinamicamente

typedef struct {

	NodoHeapCDF* nodi;              //nodi contenuti nella min-heap
	unsigned long long dimensione; //numero di nodi attualmente presenti
	unsigned long long capacita;   //numero massimo di nodi attualmente allocati

} HeapCDF;

//FUNZIONE PER INSERIRE UN RANGO NEL PATTERN COMPATTO
//rango è espresso secondo la numerazione matematica 1,...,128.
//Se rango <= 64 si modifica Erank_basso; se rango > 64 si modifica Erank_alto

void inserisci_rango_cdf(NodoHeapCDF* nodo, int rango) {

	if (rango <= 64) {

		//Il rango 1 deve corrispondere al bit 0, il rango 2 al bit 1 ecc.
		//1ULL rappresenta il valore binario 1 usando un intero unsigned long long
		//Lo shift << sposta quell'1 verso sinistra di rango-1 posizioni
		nodo->Erank_basso |= (1ULL << (rango - 1));
	}
	else {

		//Per la parte alta il rango 65 deve corrispondere al bit 0.
		//=> si sottrae 65 per ottenere l'indice del bit interno a Erank_alto
		nodo->Erank_alto |= (1ULL << (rango - 65));
	}
}


//FUNZIONE PER ELIMINARE UN RANGO DAL PATTERN COMPATTO
//Il bit corrispondente al rango indicato viene riportato a 0

void elimina_rango_cdf(NodoHeapCDF* nodo, int rango) {

	if (rango <= 64) {

		nodo->Erank_basso &= ~(1ULL << (rango - 1));
	}
	else {

		nodo->Erank_alto &= ~(1ULL << (rango - 65));
	}
}


//FUNZIONE PER VERIFICARE SE UN RANGO E' PRESENTE NEL PATTERN
//Restituisce 1 se il rango è presente e 0 se il rango non è presente

int leggi_rango_cdf(NodoHeapCDF nodo, int rango) {

	if (rango <= 64) {

		return (int)((nodo.Erank_basso >> (rango - 1)) & 1ULL);
	}
	else {

		return (int)((nodo.Erank_alto >> (rango - 65)) & 1ULL);
	}
}

//Definizione di una funzione per scambiare due nodi di una heap
//durante l'inserimento e l'estrazione di un elemento, la heap deve poter modificare la posizione dei propri nodi per mantenere in cima il pattern avente peso CDF minore.

//a e b sono puntatori ai due nodi che devono essere scambiati.

void scambia_nodi_heap(NodoHeapCDF* a, NodoHeapCDF* b) {

	NodoHeapCDF temporaneo;

	//salvo temporaneamente il primo nodo
	temporaneo = *a;

	//copio il secondo nodo nella posizione del primo
	*a = *b;

	//copio il vecchio primo nodo nella posizione del secondo
	*b = temporaneo;
}

//Definzione di una funzione per inserire un nuovo nodo nella MIN-HEAP
//La min-heap viene memorizzata tramite il vettore heap->nodi[], ma logicamente viene interpretata come un albero binario.
//Per un nodo che si trova nella posizione i:
// =>Figlio sinistro = 2*i + 1
//Figlio destro = 2*i + 2
//Padre = (i - 1) / 2

//La proprietà che deve essere rispettata è che peso del padre <= peso dei figli
//=> il nodo contenuto nella posizione 0 sarà sempre quello con il peso CDF minimo tra tutti quelli presenti


void inserisci_heap(HeapCDF* heap, NodoHeapCDF nuovo_nodo) {

	int posizione;
	int posizione_padre;

	//nel progetto N=128 la heap può contenere al massimo il numero di elementi indicato da capacita
	//Questo controllo evita di scrivere oltre i limiti della memoria allocata per nodi[]
	if (heap->dimensione >= heap->capacita) {

		printf("Errore: heap piena.\n");
		return;
	}

	//DEFINIZIONE DI UNA FUNZIONE PER ESTRARRE DALLA MIN-HEAP IL NODO AVENTE PESO CDF MINIMO
	//La proprietà della min-heap garantisce che il nodo con peso minimo si trovi sempre nella posizione 0 del vettore heap -> nodi[]
	//La funzione, quindi, salva il nodo minimo, porta l'ultimo nodo della heap in posizione 0, riduce di uno la dimensione della heap e fa scendere il nodo appena spostato fino a ripristinare la condizione necessaria => peso padre <= peso figlio
	//La funzione restituisce 1 => estrazione avvenuta correttamente oppure 0 => heap vuota, quindi non esiste alcun nodo da estrarre.
	//Il nodo estratto viene restituito attraverso il puntatore nodo_minimo.


	//il nuovo nodo viene inizialmente inserito nella prima posizione libera dell'array.
	//P.e. se dimensione=4, sono già occupate le posizioni 0,1,2,3
	//=> il nuovo nodo viene inserito nella posizione 4.
	posizione = heap->dimensione;

	heap->nodi[posizione] = nuovo_nodo;

	//dopo l'inserimento, la heap contiene un elemento in più.
	heap->dimensione++;


	//RISALITA DEL NUOVO NODO NELLA HEAP
	//Il nodo appena inserito potrebbe avere un peso minore del proprio padre. In tal caso la proprietà della min-heap sarebbe violata
	//=> il nodo viene confrontato con il padre e, se necessario, scambiato progressivamente verso l'alto

	while (posizione > 0) {

		//calcolo della posizione del padre.
		//Lavorando con indici C che partono da 0:
		//- posizione 1 o 2 => padre 0
		//- posizione 3 o 4 => padre 1
		//- posizione 5 o 6 => padre 2
		posizione_padre = (posizione - 1) / 2;

		//se il nuovo nodo ha già peso maggiore o uguale al padre, la proprietà della min-heap è rispettata
		//=> non è necessario farlo risalire ulteriormente.
		if (heap->nodi[posizione].peso_cdf >=
			heap->nodi[posizione_padre].peso_cdf) {

			break;
		}

		//se invece il nodo ha peso minore del padre, i due nodi devono essere scambiati

		scambia_nodi_heap(&heap->nodi[posizione], &heap->nodi[posizione_padre]);

		//dopo lo scambio il nuovo nodo si trova nella posizione del padre => aggiorno posizione e continuo eventualmente la risalita.

		posizione = posizione_padre;
	}
}

int estrai_minimo_heap(HeapCDF* heap, NodoHeapCDF* nodo_minimo) {

	int posizione;
	int figlio_sinistro;
	int figlio_destro;
	int posizione_minore;


	//CONTROLLO SE LA HEAP E' VUOTA
	//dimensione=0 significa che non è presente alcun nodo => non esiste un minimo da estrarre
	if (heap->dimensione == 0) {

		return 0;
	}


	//SALVATAGGIO DEL NODO MINIMO
	//per definizione della min-heap, il nodo contenuto nella posizione 0 è quello avente peso CDF minimo
	//*nodo_minimo rappresenta il nodo effettivo indicato dal puntatore passato alla funzione
	*nodo_minimo = heap->nodi[0];


	//RIDUZIONE DELLA DIMENSIONE DELLA HEAP
	//P.e. prima dell'estrazione si può avere dimensione = 5 e posizioni occupate = 0,1,2,3,4
	//dopo l'istruzione, invece, dimensione = 4 => l'ultimo nodo della vecchia heap era nella posizione 4, che ora coincide proprio con heap->nodi[heap->dimensione]

	heap->dimensione--;


	//se il nodo estratto era l'unico elemento presente, dopo la riduzione dimensione=0
	//=> la heap è diventata vuota e non occorre effettuare alcuna operazione di riordinamento
	if (heap->dimensione == 0) {

		return 1;
	}


	//SPOSTAMENTO DELL'ULTIMO NODO NELLA POSIZIONE 0
	//il minimo è stato rimosso, quindi la posizione 0 deve essere riempita.
	//Si prende l'ultimo nodo della vecchia heap e lo si porta in cima.
	//Questo nodo potrebbe però avere peso maggiore dei propri figli => sarà necessario farlo scendere
	heap->nodi[0] = heap->nodi[heap->dimensione];


	posizione = 0; //la discesa parte dalla posizione 0



	//DISCESA DEL NODO NELLA HEAP
	//ad ogni iterazione vengono individuatifigli destro, sinistro e nodo con peso minore tra padri e figli
	//Se il padre è già il più piccolo, la proprietà della heap è rispettata e la discesa termina
	while (1) {

		//per un nodo nella posizione i:
		//figlio sinistro = 2*i + 1
		//figlio destro = 2*i + 2
		figlio_sinistro = 2 * posizione + 1;
		figlio_destro = 2 * posizione + 2;


		posizione_minore = posizione; //inizialmente suppongo che il nodo corrente sia già quello avente peso minore



		//CONTROLLO DEL FIGLIO SINISTRO
		//prima verifico che il figlio sinistro esista, quindi che figlio_sinistro < heap->dimensione
		//Se esiste e ha peso minore del nodo attualmente considerato minimo, aggiorno posizione_minore

		if (figlio_sinistro < heap->dimensione &&
			heap->nodi[figlio_sinistro].peso_cdf <
			heap->nodi[posizione_minore].peso_cdf) {

			posizione_minore = figlio_sinistro;
		}


		//CONTROLLO DEL FIGLIO DESTRO
		//si effettua lo stesso controllo sul figlio destro. E', però, importante confrontarlo con posizione_minore, che potrebbe già corrispondere al figlio sinistro

		//=> alla fine posizione_minore identifica realmenteil nodo con peso più piccolo tra padre e due figli

		if (figlio_destro < heap->dimensione &&
			heap->nodi[figlio_destro].peso_cdf <
			heap->nodi[posizione_minore].peso_cdf) {

			posizione_minore = figlio_destro;
		}


		//se posizione_minore coincide ancora con posizione, significa che il nodo corrente ha peso <= a entrambi i figli
		//=> la proprietà della min-heap è rispettata e non è necessario farlo scendere ulteriormente

		if (posizione_minore == posizione) {

			break;
		}


		//se invece uno dei figli ha peso minore, il nodo corrente viene scambiato con il figlio più piccolo

		scambia_nodi_heap(&heap->nodi[posizione], &heap->nodi[posizione_minore]);


		//dopo lo scambio, il nodo che sto facendo scendere si trova nella vecchia posizione del figlio
		//=> aggiorno posizione e continuo il controllo al livello successivo

		posizione = posizione_minore;
	}


	//l'estrazione è stata completata correttamente
	return 1;
}

//Definizione di una funzione che genera i nuovi pattern CDF a partire da un nodo già estratto dalla min-heap.
//Dato un pattern con rango massimo j, vengono costruiti al massimo due nuovi pattern
//1) si aggiunge il rango j+1 mantenendo tutti i ranghi già presenti;
//2) si sostituisce il rango massimo j con il rango j+1.
//I due nuovi pattern vengono poi inseriti nella min-heap in base al rispettivo peso CDF.

void genera_nuovi_pattern_cdf(NodoHeapCDF nodo_corrente, double pesi_cdf[], HeapCDF* heap) {

	int rango_massimo;
	NodoHeapCDF pattern_aggiunta;
	NodoHeapCDF pattern_sostituzione;

	//rango_massimo è espresso secondo la numerazione matematica 1,...,N e non secondo gli indici C 0,...,N-1.
	rango_massimo = nodo_corrente.rango_massimo;

	//Se il rango massimo è già N, non esiste il rango successivo N+1 => dal nodo corrente non può essere generato alcun nuovo pattern.
	if (rango_massimo >= N) {
		return;
	}

	//I due nuovi pattern partono come copie esatte del pattern corrente.
	pattern_aggiunta = nodo_corrente;
	pattern_sostituzione = nodo_corrente;


	//PRIMO NUOVO PATTERN: AGGIUNTA DEL RANGO SUCCESSIVO
	//Se il rango massimo corrente è j, il rango j+1 viene aggiunto al pattern senza eliminare alcun rango già presente.

	inserisci_rango_cdf(&pattern_aggiunta, rango_massimo + 1);

	//Il nuovo peso si ottiene aggiungendo al peso corrente il peso CDF associato al rango j+1.
	pattern_aggiunta.peso_cdf =
		nodo_corrente.peso_cdf + pesi_cdf[rango_massimo];

	pattern_aggiunta.rango_massimo = rango_massimo + 1;

	//Il nuovo pattern viene inserito nella min-heap.
	inserisci_heap(heap, pattern_aggiunta);


	//SECONDO NUOVO PATTERN: SOSTITUZIONE DEL RANGO MASSIMO
	//Il rango massimo j viene eliminato e al suo posto viene inserito il rango j+1.

	elimina_rango_cdf(&pattern_sostituzione, rango_massimo);
	inserisci_rango_cdf(&pattern_sostituzione, rango_massimo + 1);

	//Dal peso corrente viene tolto il contributo del rango j e viene aggiunto quello del rango j+1.
	pattern_sostituzione.peso_cdf =
		nodo_corrente.peso_cdf
		- pesi_cdf[rango_massimo - 1]
		+ pesi_cdf[rango_massimo];

	pattern_sostituzione.rango_massimo = rango_massimo + 1;

	//Anche il secondo nuovo pattern viene inserito nella min-heap.
	inserisci_heap(heap, pattern_sostituzione);
}

//flag che indica se CDF-ORBGRAND ha esaurito il numero massimo di query
int abbandono = 0;

//definizione di una funzione per l'eliminazione gaussiana, prende in ingresso una matrice (A), due indici di riga (r1 e r2), senza restituire nulla (void)
//questa funzione serve a scambiare due righe della matrice
void scambia_righe(int A[NH][N], int r1, int r2) {
	int c; //c serve per scorrere le colonne
	int temp; //temp è una variabile temporanea per lo swap, che serve per non predere un valore durante lo scambio

	//nel ciclo for parto dalla colonna 0, fino a N-1, scorrendo tutte le colonne
	for (c = 0;c < N;c++) {
		temp = A[r1][c]; //salvo temporaneamente l'elemento della riga r1 e colonna c
		A[r1][c] = A[r2][c]; //metto dentro la riga r1 l'elemento che stava nella riga r2
		A[r2][c] = temp; //rimetto nella riga r2 il valore vecchio della riga r1, così le due righe si scambiano
	}
}

//riduzione generale di una matrice binaria tramite eliminazione gaussiana su GF(2)
//la funzione modifica direttamente L, memorizza in pivot_cols le colonne in cui trova i pivot
//e restituisce il numero totale di pivot trovati, cioè il rango della matrice
//È praticamente lo stesso principio già usato nella vecchia eliminazione gaussiana, solo che ora non è più vincolato alla particolare struttura del (5,2)

int riduci_gf2(int L[NH][N], int pivot_cols[NH]) {

	int r_pivot = 0;      //riga nella quale deve essere posizionato il prossimo pivot
	int c_pivot;          //colonna che si sta analizzando per cercare il pivot
	int r;                //indice di riga
	int c;                //indice di colonna
	int r_trovata;        //riga nella quale viene trovato un 1 utilizzabile come pivot

	//scorro tutte le colonne della matrice da sinistra verso destra
	//in forma_sistematica() c'era for (c_pivot = K; c_pivot < N && r_pivot < NH; c_pivot++) perché si cercavano i pivot solo nella parte destra della matrice
	//qui invece si ha for (c_pivot = 0; c_pivot < N && r_pivot < NH; c_pivot++) quindi inizio dalla prima colonna e analizzo tutte le 128 colonne
	//insomma, invece di cercare l'identità, si cercano le colonne linearmente indipendenti
	for (c_pivot = 0; c_pivot < N && r_pivot < NH; c_pivot++) {

		r_trovata = -1; //all'inizio di ogni nuova colonna
		// -1 è un valore di riferimento che funge da flag per indicare che nello specifico momento non è stata trovata nessuna riga valida

		//cerco un 1 nella colonna corrente, partendo dalla riga pivot corrente
		for (r = r_pivot; r < NH; r++) {

			if (L[r][c_pivot] == 1) { //appena trovata una riga valida, si aggiorna la label r_trovata, se alla fine è ancora -1, significa che quella colonna non può ospitare un nuovo pivot. => non si restituisce errore come nella vecchia funzione
				r_trovata = r;
				break;
			}
		}

		//se nella parte restante della colonna non è presente alcun 1,  =>questa colonna non può essere utilizzata come colonna pivot
		if (r_trovata == -1) {
			continue; //se la colonna non aggiunge un nuovo pivot => si passa alla successiva
		}

		//se il pivot è stato trovato in una riga diversa da quella desiderata, => porto la riga trovata nella posizione r_pivot
		if (r_trovata != r_pivot) {
			scambia_righe(L, r_pivot, r_trovata); //richiamo della funzione scambia_righe esattamente come in forma_sistematica()
		}

		//annullo tutti gli altri 1 presenti nella colonna pivot, sia sopra sia sotto il pivot stesso

		for (r = 0; r < NH; r++) { //scorre tutte le righe e non soltanto quelle sotto r_pivot

			if (r != r_pivot && L[r][c_pivot] == 1) {

				for (c = 0; c < N; c++) {
					L[r][c] = (L[r][c] + L[r_pivot][c]) % 2; // Per annullare un altro 1, 1+1=0 (mod-2) => si ottiene la forma desiderata
				}
			}
		}

		//memorizzo quale colonna è stata utilizzata come pivot
		//analogo a una permutazione, => per esempio, se i primi pivot si trovano nelle colonne 1, 4, 5 => pivot_cols[0]=1, pivot_cols[1]=4, pivot_col[3]=5
		//Tutte le colonne non presenti in pivot_cols saranno invece le colonne libere, da cui partità la costruzione di G
		pivot_cols[r_pivot] = c_pivot;

		//il prossimo pivot dovrà essere collocato nella riga successiva
		r_pivot++;
	}

	//il numero di righe pivot trovate coincide con il rango della matrice
	return r_pivot;
}

//definizione di una funzione che serve a portare una sequenza nella forma 000...0
//utile in corso d'opera per la sequenza Ecan, perchéogni ricerca deve partire da un pattern pulito
//V[] = vettore da azzerare; dim = la sua dimensione
void azzera_vettore(int V[], int dim) {
	int i;

	for (i = 0;i < dim;i++) { //il ciclo scorre tutte le posizioni del vettore e mette ogni elemento a zero, altrimenti Ecan potrebbe contenere valori dei tentativi precedenti
		V[i] = 0;
	}
}

//definizione di una funzione che non restituisce nulla, che serve a copiare il contenuto di un vettore in un altro
//utile se si trova una soluzione valida e si vuole salvare la code word finale in Xcan e il pattern finale in Efin
//destinazione[] = vettore dove salvare i valori; dorgente[] = vettoreda cui si prendono, dim  = numero di elementi da copiare
void copia_vettore(int destinazione[], int sorgente[], int dim) { //per ogni posizione si prende l'elemento sorgente e si salva in destinazione
	int i;

	for (i = 0;i < dim;i++) {
		destinazione[i] = sorgente[i];
	}
}

//definizione di una funzione che realizza il calcolo della sindrome SZ = ZH^T
//Z[] = parola candidata; H[NH][N] = matrice originale di controllo di parità
//SZ[NH] = vettore sindrome da riempire
void calcola_sindrome(int Z[], int H[NH][N], int SZ[NH]) {
	int r; //indice di riga
	int c; //indice di colonna

	//ogni riga di H rappresenta un vincolo di parità, quindi ogni riga produce una componente della sindrome
	for (r = 0;r < NH;r++) {
		SZ[r] = 0; //prima di accumulare i contributi della riga r-esima, si azzera a componente della sindrome, perché altrimenti dentro SZ[r] ci possono essere residui di calcolo precedente
		for (c = 0;c < N;c++) {
			SZ[r] = (SZ[r] + Z[c] * H[r][c]) % 2; //prodotto tra il simbolo di Z in posizione c per il corrispondente simbolo della riga r della matrice
		}
	}
}

//NUOVA VERSIONE SPARSA DEL CALCOLO DELLA SINDROME
//sfrutta il fatto che ogni riga della matrice LDPC H contiene esattamente 8 elementi pari a 1
void calcola_sindrome_sparsa(int Z[], int SZ_sparsa[NH]) {

	int r;
	int j;

	for (r = 0; r < NH; r++) {

		SZ_sparsa[r] = 0;

		//invece di scorrere tutte le N=128 colonne considero soltanto le 8 posizioni in cui la riga di H contiene un 1,
		for (j = 0; j < 8; j++) {

			SZ_sparsa[r] =
				(SZ_sparsa[r] + Z[H_pos[r][j]]) % 2;
		}
	}
}

//definizione di una funzione che serve a controllare se la sindrome è nulla, restituisce 0 o 1, in base al responso (1=sì, 0=no, logica booleana)
int sindrome_nulla(int SZ[NH]) {
	int r; //indice che scorre le componenti del vettore sindrome

	//scorrendo le componenti della sindrome, appena trova anche solo una componenente != 0, so già che la sindrome non è nulla => return 0
	for (r = 0;r < NH;r++) {
		if (SZ[r] != 0) {
			return 0;
		}
	}

	return 1; //se si esce dal ciclo senza aver trovato componenti non nulle, allora tutte sono 0
} //else superfluo, perché basta che il ciclo finisca senza trovare alcun 1 affinché la sindrome sia nulla e non servono condizioni

//definizione  di una funzione che prende un pattern candidato già costruito e lo mette alla prova vedendo se sia quello giusto
//Ecan[] = pattern di errore candidato; Xcan[] = code word finale salvata dal decoder
//Efin[] = pattern finale salvato dal decoder; *trovato = puntatore alla variabile che indica se il decoder ha già trovato una soluzione
void prova_pattern(int Yhard[], int H[NH][N], int Ecan[], int Z[], int SZ[NH], int Xcan[], int Efin[], int* trovato) {
	int i;

	//costruzione della candidata Z = Y ^ Ecan, per togliere al ricevuto il rumore ipotizzato, tramite XOR
	for (i = 0; i < N; i++) {
		Z[i] = Yhard[i] ^ Ecan[i];
	}


	//CONTROLLO DEL NUMERO MASSIMO DI QUERY

	//prima di effettuare una nuova query verifico se il numero massimo
	//di tentativi consentiti è già stato raggiunto
	if (numero_query >= MAX_QUERY) {

		abbandono = 1; //il decoder deve interrompere la ricerca
		return;
	}

	//NUOVO RISPETTO AL VECCHIO PROGETTO:
	//ogni volta che una candidata Z viene realmente sottoposta al controllo di sindrome, CDF-ORBGRAND ha effettuato una nuova query
	numero_query++;

	//calcolo della sindrome della candidata
	//una volta costruita Z, va testata, tramite sindrome
	calcola_sindrome_sparsa(Z, SZ);

	//se la sindrome è nulla, la candidata appartiene al codice
	if (sindrome_nulla(SZ) == 1) {
		*trovato = 1; //aggiorno il flag di successo del decoder
		copia_vettore(Xcan, Z, N); //salvo la candidata come codeword finale stimata
		copia_vettore(Efin, Ecan, N); //salvo il pattern di errore che ha portato al successo
	}
}

//Definizione di una funzione per trasformare il pattern di rango nel pattern fisico
void mappa_pattern(int Erank[], int Ecan[], int pi[]) { //Riceve il pattern dei ranghi e quello fisico da costruire insieme alla permutazione \pi, non restituendo nulla perché serve a modificare Ecan

	int k; //indice tipico del rango, k=0=>rango 1, k=1 => rango 2 [...]

	azzera_vettore(Ecan, N); //prima di costruire un nuovo pattern voglio Ecan=0...000, altrimenti si rischia di conservare gli uni appartenenti al pattern precedente

	for (k = 0;k < N;k++) { //ciclo for per scorrere i ranghi

		if (Erank[k] == 1) { //=> P.e se si ha Erank=11...000, l'if sarà vero per k=0 e k=1 e altri

			Ecan[pi[k]] = 1; //trasformazione da rango a posizione fisica
			//P.e con pi= 3 2 1... 4 0 => in realtà (4,3,2,...,5,1), per k=0 => pi[0]=3 => Ecan[3]=1=> mette 1 in posizione 3 (in C, => posizione 4 matematica)
			//Poi K=1 => pi[1]=2 => Ecan[2]=1 => Ecan = 0 0 1 1 0 
		}
	}
}

//Definizione di una funzione per trasformare il pattern CDF compatto nel pattern fisico
//La logica è la stessa di mappa_pattern(), ma il pattern dei ranghi non è più memorizzato nel vettore Erank[]
//=> la presenza del singolo rango viene verificata tramite leggi_rango_cdf()

void mappa_pattern_cdf(NodoHeapCDF nodo, int Ecan[], int pi[]) {

	int k; //indice del rango, k=0 => rango 1, k=1 => rango 2 [...]

	azzera_vettore(Ecan, N); //prima di costruire il nuovo pattern fisico si azzera Ecan

	for (k = 0; k < N; k++) {

		//leggi_rango_cdf() utilizza la numerazione matematica dei ranghi 1,...,N
		//=> al rango associato all'indice C k corrisponde k+1
		if (leggi_rango_cdf(nodo, k + 1) == 1) {

			Ecan[pi[k]] = 1; //trasformazione dal rango di affidabilità alla posizione fisica
		}
	}
}
void genera_pattern_logistici(int peso_rimanente, int rango_iniziale, int Erank[], int Ecan[], int pi[], int Yhard[], int H[NH][N], int Z[], int SZ[NH], int Xcan[], int Efin[], int* trovato) {
	

	//peso_rimanente è una variabile che indica lo stato della generazione, p.e. se all'inizio si ha peso_rimanente=5, scelto il rango 1 si avrà 5-1=4, poi, scelto il rango 4 => peso_rimanente==0, a quel punto vi sarà una combinazione valida
	//rango_iniziale --> serve a impedire di utilizzare due volte lo stesso rango ed evitare duplicati
	int rango;

	//if (*trovato == 1) {
		//return;
	//}

	//MODIFICA: la ricorsione deve arrestarsi sia quando viene trovata una codeword valida sia quando viene raggiunto il numero massimo di query
	if (*trovato == 1 || abbandono == 1) {
		return;
	}

	if (peso_rimanente == 0) { //il logistic Weight sarà stato costruito completamente, Erank contiene quindi un pattern valido nel dominio dei ranghi, => lo trasformo nelle corrispondenti posizioni fisiche tramite \pi

		mappa_pattern(Erank, Ecan, pi); //=> attuo la trasformazione dal dominio dei ranghi alle posizioni fisiche tramite \pi

		prova_pattern(Yhard, H, Ecan, Z, SZ, Xcan, Efin, trovato);

		return;
	}

	for (rango = rango_iniziale; rango <= N; rango++) { //Ciclo for con cui attuo uno alla volta tutti i ranghi ancora utilizzabili
		//se scelgo un rango, lo inserisco nel pattern, sottraggo il suo valore dal Logistic Weight ancora da costruire, provo ricorsivamente a completare il pattern utilizzando soltanto ranghi successivi e, terminato quel tentativo, annullo la scelta per provarne un'altra
		//qui rango rappresenta il rango di affidabilità che stiamo valutando come possibile elemento del pattern. Non parte sempre da 1, ma da rango_iniziale

		//Un rango può essere selezionato solo se il suo contributo al Logistic Weight non supera il peso ancora da costruire. Se, p.e., peso_rimanente=3, non avrebbe senso scegliere rango=4, perché la somma supererebbe già il Logistic Weight desiderato.

		if (rango <= peso_rimanente) {

			Erank[rango - 1] = 1; //Seleziono provvisoriamente il rango corrente inserendo un 1 nella corrispondente posizione di Erank. Poiché i ranghi matematici vanno da 1 a N mentre gli indici C vanno da 0 a N-1, il rango r viene memorizzato nella posizione r-1.

			//Una volta scelto il rango corrente, richiamo ricorsivamente la funzione per completare il pattern. Il peso ancora da costruire viene ridotto del valore del rango scelto: nuovo_peso = peso_rimanente - rango.
			//La ricerca successiva parte da rango + 1, così lo stesso rango non può essere utilizzato nuovamente e le combinazioni vengono costruite in ordine strettamente crescente, evitando duplicati.
			//Erank mantiene memoria dei ranghi già selezionati durante il percorso ricorsivo; Ecan e pi serviranno per la successiva conversione del
			//pattern di rango nelle effettive posizioni fisiche.

			genera_pattern_logistici(peso_rimanente - rango, rango + 1, Erank, Ecan, pi, Yhard, H, Z, SZ, Xcan, Efin, trovato);

			Erank[rango - 1] = 0; //Terminata l'esplorazione ricorsiva di tutti i pattern che contengono il rango corrente, annullo la scelta riportando a 0 la corrispondente posizione di Erank
			//Questo è il passo di backtracking: consente di ritornare allo stato precedente del pattern e provare, nell'iterazione successiva del for,
			//un rango differente senza conservare gli 1 appartenenti al ramo appena esplorato

			//if (*trovato == 1) {
				//return;
			//} //la stessa funzione in "AWGN" si trasforma da generatore di prova a generatore del decoder vero

			//MODIFICA: arresto del backtracking anche in caso di abbandono
			if (*trovato == 1 || abbandono == 1) {
				return;
			}
		}
	}
}

//Definizione della CDF della distribuzione gaussiana standard
//erf() è la funzione errore, utilizzabile grazie alla libreria math.h, che permette di calcolare Phi(x) = P(Z<=x) con Z distribuita ovviamente da N(0,1)
double cdf_normale_standard(double x) {
	return 0.5 * (1.0 + erf(x / sqrt(2.0)));
}

//Definizione della CDF Psi(t) del modulo degli LLR per il canale AWGN
//t=possibile valore del modulo di LLR
double cdf_modulo_llr(double t, double Es, double sigma) {

	double mu_llr;
	double sigma_llr;

	//l' LLR condizionato a uno dei due simboli della BPSK sarà gaussiano con media = 2Es/sigma quadro

	mu_llr = (2.0 * Es) / (sigma * sigma);

	//la deviazione standard, invece, = 2\sqrt(es)/sigma
	sigma_llr = (2.0 * sqrt(Es)) / sigma;

	//Psi(t)=P(|LLR|<=t) => La seconda Phi deriva dal limite inferiore -t del valore assoluto.
	return cdf_normale_standard((t - mu_llr) / sigma_llr) + cdf_normale_standard((t + mu_llr) / sigma_llr) - 1.0;
}

//Definizione dell'inversa numerica della CDF del modulo dell'LLR.
//Data una probabilità compresa tra 0 e 1, la funzione ricerca il valore t tale che cdf_modulo_llr(t)=probabilita.
//La ricerca viene effettuata tramite bisezione, restringendo progressivamente un intervallo che contiene la soluzione.

double inversa_cdf_modulo_llr(double probabilita, double Es, double sigma) {

	double limite_inferiore = 0.0;
	double limite_superiore = 100.0;
	double valore_medio;
	double probabilita_media;
	int iterazione;

	//Si eseguono 100 iterazioni, numero più che sufficiente per ottenere una precisione elevata nel calcolo del peso CDF.
	for (iterazione = 0; iterazione < 100; iterazione++) {

		valore_medio = (limite_inferiore + limite_superiore) / 2.0;

		probabilita_media = cdf_modulo_llr(valore_medio, Es, sigma);

		//Se la probabilità ottenuta è ancora minore di quella desiderata, il valore cercato deve trovarsi nella metà superiore dell'intervallo.
		if (probabilita_media < probabilita) {
			limite_inferiore = valore_medio;
		}
		else {
			limite_superiore = valore_medio;
		}
	}

	//Al termine della bisezione, i due estremi risultano praticamente coincidenti; il loro punto medio rappresenta l'approssimazione dell'inversa della CDF.
	return (limite_inferiore + limite_superiore) / 2.0;
}
//Definisco una funzione che rstituisce un valore reale. In questa fase voglio ottenere una normale standard, => media nulla e varianza 1.
double genera_gaussiana_standard(void) {

	const double PI = 3.14159265358979323846;


	double U1; //generazione degli ingressi della funzione
	double U2;
	double Z;

	//generazione di due variabili uniformi nell'intervallo (0,1)

	U1 = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0); //la presenza di +1.0 e +2.0 sono importanti per ottenere calori strettamente compresi tra 0 e 1
	U2 = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);

	//in particolare è importante per U1, infatti ln(U1) fa si che U1 debba essere !=0

	//trasformazione di Box-Muller, che trasforma due variabili uniformi in una normale

	Z = sqrt(-2.0 * log(U1)) * cos(2.0 * PI * U2);

	return Z; //per restituire il campione al main
}

//definizione di una funzione per ordinare le affidabilità, universale per poter essere usata anche a progetto completo
void ordina_affidabilita(double Affidabilita[], int pi[]) { //riceve i due vettori Affidabilita[] e pi[], è void dunque modifica i valori di pi[]


	// Affidabilita[] contiene i moduli degli LLR e verrà soltsnto letto, mentre pi[] gli indici delle posizioni fisiche e verrà modificato
	int i; // Posizione di pi che va sistemata
	int j; //serve per cercare, nelle posizioni successive, l'affidabilità più piccola
	int indice_min; //Memorizza la posizione dentro pi nella quale, fino a quel momento, si ha trovato l'affidabilità minima
	int temp; //per scambiare elementi

	for (i = 0;i < N - 1;i++) { //=> selection sort, va decisio quale indice mettere in pi[i], si cerca tra p[i], pi[i+1] ecc..., quello associato all'affidabilità più piccola

		// si parte da i=0, => va trovato il bit meno affidabile di tutti e metterne l'indice in pi[0]
		//poi per i=1 si cercherà il meno affidabile tra quelli restanti e così via
		// //c'è N-1 perché l'ultimo bit rimasto non ha bisogno di cercare posizione

		indice_min = i; //suppongo che il minimo sia quello che si trova nella poszione i di pi. => poi controllo gli altri e vedo se ne trovo uno più piccolo

		//si controllano tutti gli elementi successivi
		//Dato che si ha già assunto quello in poszione i come minimo provvisorio, si parte da J=i+1
		// con i=0 si ha j=1,2,3,4 => va confrontato l'elemento attuale con ciò che segue
		for (j = i + 1;j < N;j++) {
			if (Affidabilita[pi[j]] < Affidabilita[pi[indice_min]]) { // non si confronta l'affidabilità, perché pi potrebbe essere già stato riordinato.

				indice_min = j; //si seguono gli indici contenuti in pi
				//P.e pi = 0 1 2 3 4 ,  per j=1 si ha pi[j]=pi[1]=1; => Affidabilita[pi[j]] diventa Affidabilita[1]
				//Dall'altra parte indice_min=0 => pi[indice_min]=pi[0]=0 e Affidabilita[pi[indice_min]] diventa Affidabilita[0]
				// nell'esempio si aveva Affidabilità[pi[j]] => Affidabilita[1] = 1.90 e Affidabilita[pi[indice_min]] => Affidabilita[0]=4.60
				// => Chiedo 1.90 < 4.60 ? sì, allora eseguo indice_min=j => indice_min =1 => il minimo trovato non è più in poszione 0 di pi, ma 1
				//itero il processo
			}
		}

		//qui chiedo se il valore minimo trovato era già nella posizione corretta o va spostato
		//P.e se si aveva i=0 e poi indice_min=4 => 4 != 0 => va spostato, lo scambio si fa con la variabile temp
		if (indice_min != i) {

			temp = pi[i];
			pi[i] = pi[indice_min];
			pi[indice_min] = temp;
		}
	}
	//qui si usa una semplice selection sort, ma non si confronta Affidabilità[j] direttamente
	//si confronta Affidabilità[pi[j]], perché si riordinano gli indici, non il vettore di affidabilità
	//in altre parole non voglio modificare il vettore delle affidabilità; voglio ordinare gli indici delle posizioni fisiche dei bit in base alle rispettive affidabilità
	//occorre, d'altra parte, fare attenzione agli indici in C; infatti il C utilizza lo 0, non si avrà p.e 1,2,3,4,5 ; ma 0,1,2,3,4 .
	// => nella stampa finale si aggiungerà pi[i] + 1;
}


int main() {

	int H[NH][N]; //crea la matrice H (64 x 128)
	int L[NH][N]; //copia di lavoro di H, da modificare durante l'eliminazione gaussiana

	int pivot_cols[NH]; //memorizza le colonne in cui vengono individuati i pivot
	int rango;          //numero di pivot trovati = rango della matrice H

	int free_cols[N]; //conterrà gli indici delle colonne non pivot
	int num_free;     //numero effettivo di colonne libere
	int K;            //dimensione del codice

	int r;
	int c;

	FILE* fp = fopen("H_LDPC_128_wc4_wr8.txt", "r"); //apertura del file di testo estratto da MATLAB; "r" indica la modalità lettura (read)

	if (fp == NULL) {
		printf("Errore nell'apertura del file contenente H.\n");
		return 1;
	}

	//doppio ciclo che visita tutte le celle della matrice => 64*128=8192
	for (r = 0; r < NH; r++) {

		for (c = 0; c < N; c++) {

			//fscan deve riuscire a leggere un intero per ogni posizione; => se non ci riesce, il programma dice esattamente dove si è fermato tramite if
			if (fscanf(fp, "%d", &H[r][c]) != 1) {

				printf("Errore nella lettura di H alla riga %d, colonna %d.\n", r, c);

				fclose(fp);
				return 1;
			}
		}
	}

	fclose(fp);

	//copio H nella matrice di lavoro L
	for (r = 0; r < NH; r++) {

		for (c = 0; c < N; c++) {

			L[r][c] = H[r][c];
		}
	}

	printf("Matrice H letta correttamente.\n");
	printf("Dimensioni previste: %d x %d\n", NH, N);

	int peso_riga;

	for (r = 0; r < NH; r++) {

		peso_riga = 0; //azzero il contatore

		for (c = 0; c < N; c++) {
			peso_riga = peso_riga + H[r][c]; //somma di tutti i valori della riga; dato che H è binaria, questa somma equivale esattamente al numero di 1 presenti nella riga, cioè al suo peso di Hamming

			//la costruzione di Gallager usata impone come peso di riga 8, => per ognuna delle 64 righe, vanno trovati precisamente 8 uni.
		}

		if (peso_riga != 8) { //=> se anche solo una riga non soddisfa questo vincolo, il programma segnala l'anomalia
			printf("Errore: la riga %d ha peso %d invece di 8.\n", r, peso_riga);
			return 1;
		}
	}

	printf("Controllo righe superato: tutte le righe hanno peso 8.\n");

	//COSTRUZIONE DELLA RAPPRESENTAZIONE SPARSA DI H, serve principalmente per ottimizzare la resa della compilazione dal punto di vista del tempo di elaborazione

	for (r = 0; r < NH; r++) {

		int indice_uno = 0;

		for (c = 0; c < N; c++) {

			if (H[r][c] == 1) {

				H_pos[r][indice_uno] = c;
				indice_uno++;
			}
		}
	}

	printf("Rappresentazione sparsa di H costruita correttamente.\n");

	int peso_colonna;

	for (c = 0; c < N; c++) {

		peso_colonna = 0; //azzero il contatore per ogni nuova colonna

		for (r = 0; r < NH; r++) {
			peso_colonna = peso_colonna + H[r][c];
			//sommo tutti gli elementi della colonna;
			//essendo H binaria, la somma coincide con il numero di 1 presenti
		}

		//la costruzione di Gallager impone peso di colonna pari a 4
		if (peso_colonna != 4) {
			printf("Errore: la colonna %d ha peso %d invece di 4.\n",
				c, peso_colonna);
			return 1;
		}
	}

	printf("Controllo colonne superato: tutte le colonne hanno peso 4.\n");

	//riduzione della matrice di lavoro L su GF(2)
	//la funzione restituisce il numero di pivot trovati, cioè il rango di H
	rango = riduci_gf2(L, pivot_cols);

	printf("Rango di H calcolato in C: %d\n", rango);

	K = N - rango; //istruzione immediata, implementa il calcolo della lunghezza del codice

	printf("Dimensione del codice K calcolata in C: %d\n", K);

	int is_pivot[N] = { 0 }; //questa inizializzazione porta tutti i 128 elementi a zero.
	//=> all'inizio considero tutte le colonne non pivot
	//Poi marco quelle che riduci_gf2() ha effettivamente trovato


	//Marco con 1 tutte le colonne che riduci_gf2() ha individuato come pivot
	//quando r = 0, prende pivot-cols[0] che contiene l'indice della prima colonna pivot, e si fa is_pivot[pivot_cols[0]]=1
	//Poi per r=1 lo stesso, is_pivot[pivot_cols[1]]=1 eccetera. Il ciclo arriva fino a r<rango => 61 valori
	for (r = 0; r < rango; r++) {
		is_pivot[pivot_cols[r]] = 1;
	}

	//Cerco le colonne che non sono pivot e ne salvo gli indici in free_cols
	//=> operazione opposta
	//Contatore di colonne libere trovate
	num_free = 0;

	for (c = 0; c < N; c++) {

		if (is_pivot[c] == 0) {

			free_cols[num_free] = c;
			num_free++;
		}
	}

	printf("Numero di colonne pivot: %d\n", rango);
	printf("Numero di colonne libere: %d\n", num_free);

	//il numero di colonne libere deve coincidere con la dimensione K del codice
	if (num_free != K) {
		printf("Errore: il numero di colonne libere (%d) non coincide con K (%d).\n", num_free, K);
		return 1;
	}

	printf("Controllo superato: il numero di colonne libere coincide con K.\n");

	//allocazione dinamica della matrice generatrice G
	//G avrà K righe e N colonne, quindi nel caso corrente dimensione 67 x 128
	//G è un puantatore a righe costituite ciascuna da N interi, poiché N=128, ogni riga puntata da G contiene 128 interi 
	int (*G)[N] = calloc(K, sizeof * G); //calloc riserva memoria sufficiente per K righe della stessa dimensionedi una riga di G => verrà allocata una matrice 67x128

	if (G == NULL) {
		printf("Errore nell'allocazione della matrice generatrice G.\n");
		return 1;
	}

	int q;//indice della riga di G, associata a una variabile libera
	int p;//indice utilizzato per scorrere le colonne pivot
	int free_col;//indice fisico della colonna libera corrente

	//costruzione della matrice generatrice G
	//ogni riga di G viene associata ad una delle K colonne libere
	for (q = 0; q < K; q++) {

		//indice fisico della q-esima colonna libera
		free_col = free_cols[q];

		//La variabile libera associata alla riga corrente viene posta uguale a 1
		//le altre colonne libere restano a 0 grazie all'inizializzazione effettuata da calloc
		G[q][free_col] = 1;

		//Determino il valore delle variabili pivot utilizzando le equazioni contenute nella matrice ridotta L

		for (p = 0; p < rango; p++) {

			//Punto fondamentale, una volta scelta la colonna libera associata alla riga corrente di G, leggo nella matrice ridotta L quali valori devono assumere le colonne pivot e li copio nelle corrispondenti posizioni della riga di G.
			//Per costrire ogni riga della matrice generatrice considero una variabile libera alla volta: la pongo a 1, mentre tutte le altre variabili libere rimangono a 0. La matrice ridotta L consente poi di ricavare i valori delle variabili pivot necessari a soddisfare le equazioni di parità. 
			//Per questo, per ogni pivot, il valore da inserire nella riga di G coincide con l'elemento di L posto sulla stessa riga pivot e sulla colonna libera considerata.
			G[q][pivot_cols[p]] = L[p][free_col];
		}
	}

	//CONTROLLO DELLA PARTE DI G CORRISPONDENTE ALLE COLONNE LIBERE

	//durante la costruzione di G, per ogni riga q, è stata posta a 1 solamente la q-esima variabile libera
	//tutte le altre variabili libere sono invece rimaste a 0 grazie all'inizializzazione effettuata con calloc
	//=> se considero solamente le colonne contenute in free_cols, devo quindi ottenere una matrice identità KxK

	int j; //indice che serve per scorrere le K colonne libere
	int controllo_identita = 1; //Flag di controllo; parto assumendo che il controllo sia corretto
	//se durante la verifica trovo anche un solo elemento errato, questo flag viene posto a 0

	for (q = 0; q < K; q++) { //scorro tutte le K righe della matrice G

		for (j = 0; j < K; j++) { //Per ogni riga controllo tutte le K colonne libere

			//free_cols[j] contiene l'indice fisico della j-esima colonna libera nella matrice G. Nella riga q, la colonna libera associata proprio a q deve contenere 1 mentre tutte le altre colonne libere devono contenere 0
			//l'espressione (q == j) in C restituisce: 1 se q e j coincidono o 0 se q e j sono diversi
			//=> in questo modo posso confrontare direttamente l'elemento di Gcon il valore che dovrebbe assumere in una matrice identità

			if (G[q][free_cols[j]] != (q == j)) {

				//se trovo anche solo un elemento diverso da quello previsto,
				//la parte di G relativa alle colonne libere non coincide con I_K
				controllo_identita = 0;

				break; //non serve continuare a controllare le altre colonne della stessa riga
			}
		}

		//se il controllo è già fallito in una riga precedente, posso uscire anche dal ciclo esterno
		if (controllo_identita == 0) {
			break;
		}
	}

	//verifica finale del risultato del controllo
	if (controllo_identita == 1) {

		printf("Controllo G superato: sulle colonne libere e' presente I_K.\n");
	}
	else {

		printf("Errore nella costruzione di G: il blocco sulle colonne libere non e' I_K.\n");

		//G è stata allocata dinamicamente, quindi prima di terminare il programma
		//devo comunque liberare la memoria precedentemente riservata
		free(G);

		return 1;
	}

	//CONTROLLO DI ORTOGONALITA' TRA H E G

	//Perché G sia realmente una matrice generatrice del codice definito da H, deve valere H * G^T = 0 su GF(2)
	//=> ogni elemento del prodotto H*G^T deriva dal prodotto scalare tra una riga di H e una riga di G
	//=> devo verificare tutte le NH righe di H con tutte le K righe di G

	int controllo_ortogonalita = 1; //Flag; inizialmente assumo che il controllo sia corretto
	int somma; //serve per calcolare il prodotto scalare mod-2 tra una riga di H e una riga di G

	for (r = 0; r < NH; r++) { //scorro le 64 righe di H

		for (q = 0; q < K; q++) { //Per ogni riga di H considero tutte le 67 righe di G

			somma = 0; //azzero il risultato prima di ogni nuovo prodotto scalare

			for (c = 0; c < N; c++) { //scorro tutte le 128 posizioni delle due righe

				//Moltiplico gli elementi corrispondenti di H e G e sommo i contributi
				//l'operazione viene sempre mantenuta in GF(2) tramite modulo 2
				somma = (somma + H[r][c] * G[q][c]) % 2;
			}

			//se il prodotto scalare vale 1, la riga q di G non soddisfa
			//almeno uno dei vincoli di parita' imposti dalla riga r di H
			if (somma != 0) {

				controllo_ortogonalita = 0;

				printf("Errore di ortogonalita' tra la riga %d di H e la riga %d di G.\n",
					r, q);

				break;
			}
		}

		//se è già stato trovato un errore non serve continuare con le altre righe
		if (controllo_ortogonalita == 0) {
			break;
		}
	}

	//Controllo finale
	if (controllo_ortogonalita == 1) {

		printf("Controllo di ortogonalita' superato: H * G^T = 0.\n");
	}
	else {

		printf("Errore: la matrice G costruita non soddisfa H * G^T = 0.\n");

		free(G);
		return 1;
	}

	//allocazione dinamica del messaggio informativo M, formato da K bit, => si chiede memoria per K=67 interi
	int* M = malloc(K * sizeof * M);

	if (M == NULL) {
		printf("Errore nell'allocazione del messaggio informativo M.\n");
		free(G);
		return 1;
	}
	int X[N]; //sequenza codificata di lunghezza N=128

	//si inizializza il generatore pseudo-casuale
	//srand((unsigned int)time(NULL));

	//SEME PSEUDOCASUALE FISSATO PER IL CONFRONTO TRA ORBGRAND E CDF-ORBGRAND
	//Lo stesso seme viene utilizzato nei due progetti => a parità di sequenza di chiamate a rand(), vengono generati gli stessi messaggi informativi e gli stessi campioni di rumore
	//In questo modo il confronto tra i decoder viene effettuato sulle stesse realizzazioni del canale e le differenze osservate dipendono dall'ordine di ricerca dei pattern e non da differenti sequenze casuali

	srand(12345);

	//CONTATORI PER IL TEST SU PIU' BLOCCHI A Eb/N0 = 4 dB

	int numero_blocchi = 0;       //numero totale di blocchi trasmessi
	int numero_corretti = 0;      //blocchi decodificati correttamente
	int numero_undetection = 0;   //codeword valida trovata, ma diversa da X
	int numero_abbandoni = 0;     //decoder fermato per MAX_QUERY

	//numero di errori totali che si vogliono raccogliere per ogni valore di Eb/N0
	//la simulazione non viene più arrestata dopo un numero fissato di blocchi, ma continua finché la somma tra undetection e abbandoni raggiunge tale valore
	//=> TARGET_ERRORI rappresenta il criterio di arresto della simulazione
	const int TARGET_ERRORI = 100;

	unsigned long long query_totali = 0; //somma delle query effettuate su tutti i blocchi
	unsigned long long query_massime_successo = 0; //massimo numero di query osservato nei blocchi in cui CDF-ORBGRAND trova una codeword

	//CONTATORI DEDICATI A ORBGRAND NEL PROGETTO DI CONFRONTO
	//Vengono mantenuti separati dai contatori già utilizzati da CDF-ORBGRAND => i due decoder possono essere confrontati sugli stessi blocchi senza mescolare i rispettivi risultati

	int numero_corretti_orb = 0;      //blocchi decodificati correttamente da ORBGRAND
	int numero_undetection_orb = 0;   //codeword valida trovata da ORBGRAND, ma diversa da X
	int numero_abbandoni_orb = 0;     //ORBGRAND fermato per MAX_QUERY

	unsigned long long query_totali_orb = 0; //somma delle query effettuate da ORBGRAND su tutti i blocchi
	unsigned long long query_massime_successo_orb = 0; //massimo numero di query osservato nei blocchi in cui ORBGRAND trova una codeword

	//CONTATORI DEI RECUPERI NEL CONFRONTO DIRETTO TRA I DUE DECODER
	//Un recupero indica che un decoder ricostruisce correttamente X mentre l'altro raggiunge MAX_QUERY sullo stesso blocco

	int recuperi_cdf = 0; //CDF-ORBGRAND corregge un blocco abbandonato da ORBGRAND
	int recuperi_orb = 0; //ORBGRAND corregge un blocco abbandonato da CDF-ORBGRAND

	//CONTATORI DEL CONFRONTO DELLE QUERY NEI BLOCCHI CORRETTI DA ENTRAMBI I DECODER
	//Il confronto viene effettuato soltanto quando sia ORBGRAND sia CDF-ORBGRAND ricostruiscono correttamente X

	int cdf_meno_query = 0; //numero di blocchi in cui CDF-ORBGRAND corregge usando meno query di ORBGRAND
	int orb_meno_query = 0; //numero di blocchi in cui ORBGRAND corregge usando meno query di CDF-ORBGRAND
	int query_uguali = 0;   //numero di blocchi in cui entrambi correggono utilizzando lo stesso numero di query

	int i; //indice utilizzato per scorrere i K bit del messaggio informativo

	//Parametri del codice
	double Rc = (double)K / N; // double Rc = (double)K / N; //il (double) è essenziale perché K e N sono valori interi; senza il cast 67/128 verrebbe eseguita come divisione intera e darebbe 0

	//energia del simbolo BPSK normalizzata
	double Es = 1.0;

	//Rapporto Eb/N0 scelto per la simulazione
	double EbN0_dB = 7.0; //si aggiunge Eb/N0 che conduce a \sigma per analizzare il rumore del canale
	//tale per cui la variabile Z distribuita con PDF normale standard sarà esattamente un elemento con cui si definisce il rumore additivo, con n=\sigma*Z
	//Fissato Eb/N0 a 4 dB come valore di prova
	//Fissato Eb/N0 a 4 dB come valore iniziale di prova;
	//in questo caso il rate utilizzato non è più quello del codice toy,
	//Ma Rc = K/N = 67/128, ricavato dalla dimensione effettiva del codice LDPC

	//conversione da dB a unità assoluta
	double gamma_b = pow(10.0, EbN0_dB / 10.0);

	//varianza teorica del rumore AWGN
	double sigma_quadro = Es / (2.0 * Rc * gamma_b);

	//deviazione standard del rumore
	double sigma = sqrt(sigma_quadro);

	double pesi_cdf[N]; //vettore che conterrà i pesi CDF calcolati automaticamente in funzione del canale

	//CALCOLO AUTOMATICO DEI PESI CDF
	//Per ciascun rango r=1,...,N si considera la probabilità r/(N+1) e si calcola il corrispondente quantile della distribuzione del modulo dell'LLR
	//La posizione i del vettore corrisponde al rango matematico r=i+1, poiché gli indici C partono da 0

	for (i = 0; i < N; i++) {

		double probabilita_cdf = (double)(i + 1) / (N + 1);

		pesi_cdf[i] = inversa_cdf_modulo_llr(probabilita_cdf, Es, sigma);
	}

	//stampa di controllo dei pesi CDF
	printf("\n=====================================\n");
	printf("PESI CDF CALCOLATI DAL CANALE\n");
	printf("=====================================\n");

	for (i = 0; i < N; i++) {
		printf("Rango %d => peso CDF = %.6f\n", i + 1, pesi_cdf[i]);
	}

	printf("\n");

	//MIN-HEAP UTILIZZATA DAL DECODER CDF-ORBGRAND
	//La memoria viene allocata una sola volta prima della simulazione dei blocchi e viene poi riutilizzata per ogni nuova decodifica

	HeapCDF heap_decoder_cdf;

	heap_decoder_cdf.dimensione = 0;
	heap_decoder_cdf.capacita = MAX_QUERY;

	//Ogni elemento della heap è un NodoHeapCDF; viene riservato spazio sufficiente per il numero massimo di nodi previsto
	heap_decoder_cdf.nodi = malloc(heap_decoder_cdf.capacita * sizeof * heap_decoder_cdf.nodi);

	if (heap_decoder_cdf.nodi == NULL) {

		printf("Errore nell'allocazione della min-heap CDF.\n");

		free(M);
		free(G);

		return 1;
	}


	//SIMULAZIONE SU UN NUMERO VARIABILE DI BLOCCHI
	//non si simulano più 100 blocchi fissati a priori;
	//la simulazione deve continuare finché vengono osservati 100 errori totali, dove per errore totale si intende undetection + abbandono
	//=> il numero di blocchi trasmessi non è più fissato a priori, ma sarà un valore M dipendente dalle prestazioni del decoder al valore di Eb/N0 considerato


	while (numero_undetection + numero_abbandoni < TARGET_ERRORI || numero_undetection_orb + numero_abbandoni_orb < TARGET_ERRORI) {
		//ogni iterazione del while corrisponde alla trasmissione
		//e successiva decodifica di un nuovo blocco
		numero_blocchi++;

		printf("\n\n=========================================\n");
		printf("BLOCCO %d\n", numero_blocchi);
		printf("=========================================\n");

		//GENERAZIONE CASUALE DI M

		//generazione casuale dei K bit del messaggio informativo
		//ogni posizione può assumere solamente 0 oppure 1
		//=> rand()%2 restituisce direttamente un bit binario
		for (c = 0; c < K; c++) {
			M[c] = rand() % 2;
		}

		//CODIFICA TRAMITE G
		//Per costruire ogni bit della codeword scorro tutti i bit del messaggio informativo e considero il contributo della corrispondente riga di \(G\); la somma viene effettuata modulo 2.

		//CODIFICA DEL MESSAGGIO INFORMATIVO
		//calcolo X = M * G in GF(2)
		for (c = 0; c < N; c++) {

			X[c] = 0; //azzero il bit della codeword prima di calcolarlo

			for (i = 0; i < K; i++) {

				X[c] = (X[c] + M[i] * G[i][c]) % 2;
			}
		}


		if (STAMPE_DETTAGLIATE) {

			//STAMPA DI M E DI X
			printf("\nMessaggio informativo M:\n");

			for (i = 0; i < K; i++) {
				printf("%d", M[i]);
			}

			printf("\n");

			printf("\nCodeword X:\n");

			for (c = 0; c < N; c++) {
				printf("%d", X[c]);
			}

			printf("\n");
		}

		//CONTROLLO DELLA CODEWORD X TRAMITE SINDROME
		//Dopo aver codificato M tramite G, X deve essere una parola valida del codice
		//=> Deve soddisfare tutti i vincoli di parità imposti da H
		//=> calcolo la sindrome S_X = X * H^T su GF(2), che deve risultare completamente nulla

		int SX[NH]; //sindrome della codeword X, composta da NH=64 elementi
		int controllo_sindrome_X = 1; //flag; inizialmente assumo che la sindrome sia nulla

		for (r = 0; r < NH; r++) { //ogni riga di H corrisponde ad un vincolo di parità

			SX[r] = 0; //azzero l'elemento corrente della sindrome

			for (c = 0; c < N; c++) { //scorro tutti i 128 bit della codeword

				//calcolo il prodotto scalare mod-2 tra X e la riga r di H
				//=> se il risultato è 0, X soddisfa quello specifico vincolo di parità
				SX[r] = (SX[r] + X[c] * H[r][c]) % 2;
			}

			//se anche un solo elemento della sindrome vale 1,
			//la sequenza X ottenuta dalla codifica non è una parola valida del codice
			if (SX[r] != 0) {

				controllo_sindrome_X = 0;

				printf("Errore: sindrome non nulla alla posizione %d.\n", r);

				break;
			}
		}

		//verifica finale del controllo
		if (controllo_sindrome_X == 1) {

			if (STAMPE_DETTAGLIATE) {
				printf("Controllo codifica superato: la codeword X ha sindrome nulla.\n");
			}
		}
		else {

			printf("Errore nella codifica: X non soddisfa tutti i vincoli di parita' di H.\n");

			free(M);
			free(G);

			return 1;
		}

		//CONTROLLO DELLA SISTEMATICITA' DELLA CODIFICA
		//la matrice G costruita contiene una matrice identità nelle colonne libere
		//=> ogni bit informativo M[q] deve essere riportato direttamente
		//nella posizione della codeword indicata da free_cols[q]

		int controllo_sistematicita = 1; //flag; inizialmente assumo che il controllo sia corretto

		for (q = 0; q < K; q++) {

			//confronto il q-esimo bit informativo con la posizione della codeword
			//associata alla q-esima colonna libera
			if (X[free_cols[q]] != M[q]) {

				controllo_sistematicita = 0;

				printf("Errore di sistematicita': M[%d] diverso da X[%d].\n",
					q, free_cols[q]);

				break;
			}
		}

		if (controllo_sistematicita == 1) {

			if (STAMPE_DETTAGLIATE) {
				printf("Controllo sistematicita' superato: M coincide con X nelle colonne libere.\n");
			}
		}

		else {

			free(M);
			free(G);

			return 1;
		}

		//MODULAZIONE BPSK
		//la codeword X è ancora costituita da bit 0 e 1;
		//Prima della trasmissione sul canale AWGN ogni bit viene trasformato in un simbolo reale utilizzo la stessa convenzione adottata nel precedente programma: 0 -> +1 e 1 -> -1

		double Sbpsk[N]; //sequenza dei simboli BPSK corrispondenti ai N=128 bit della codeword

		for (c = 0; c < N; c++) {

			if (X[c] == 0) {
				Sbpsk[c] = 1.0;
			}
			else {
				Sbpsk[c] = -1.0;
			}
		}

		if (STAMPE_DETTAGLIATE) {

			printf("\nSequenza modulata BPSK:\n");

			for (c = 0; c < N; c++) {
				printf("%.0f ", Sbpsk[c]);
			}

			printf("\n");
		}

		printf("\n");

		//VETTORI NECESSARI PER LA TRASMISSIONE SUL CANALE AWGN
		double RumoreCanale[N]; //campioni di rumore additivo applicati ai simboli BPSK
		double Ysoft[N];        //valori ricevuti dal canale prima della decisione hard

		int Yhard[N];           //sequenza ottenuta dopo la decisione hard
		double LLR[N];          //Log-Likelihood Ratio associato ad ogni bit ricevuto
		double Affidabilita[N]; //modulo degli LLR, utilizzato per valutare l'affidabilita' dei bit
		int pi[N];              //permutazione degli indici fisici ordinati dal bit meno affidabile al piu' affidabile

		double Zgauss;          //singolo campione gaussiano standard N(0,1)

		//VARIABILI DI LAVORO UTILIZZATE DAL DECODER

		int trovato = 0; //flag che indica se il decoder ha già trovato una soluzione, all'inizio non ne è stata trovata una

		int wL; //Logistic Weight corrente

		int wL_max = N * (N + 1) / 2; //massimo Logistic Weight possibile
		//queste due variabili servono perché il decoder esplorerà automaticamente
		//wL = 0, 1, ..., wL_max

		//definizione dei vettori di lavoro utilizzati dal decoder

		int Erank[N]; //rappresenta il pattern nel dominio dei ranghi

		int Ecan[N]; //pattern candidato riportato nelle posizioni fisiche della codeword

		int Z[N]; //parola candidata ottenuta come Z = Yhard XOR Ecan

		int SZ[NH]; //sindrome associata alla parola candidata Z

		//VETTORI UTILIZZATI PER CONFRONTARE IL CALCOLO DENSO E IL CALCOLO SPARSO DELLA SINDROME
		int SZ_densa_test[NH];
		int SZ_sparsa_test[NH];

		int Xcan[N]; //codeword finale stimata dal decoder

		int Efin[N]; //pattern di errore che ha prodotto la codeword valida trovata

		//VARIABILI PER SALVARE IL RISULTATO DI ORBGRAND SUL SINGOLO BLOCCO
		//Sono necessarie perché, terminato ORBGRAND, CDF-ORBGRAND riutilizzerà e azzererà numero_query, trovato e abbandono

		unsigned long long query_blocco_orb = 0; //numero di query effettuate da ORBGRAND sul blocco corrente

		int corretto_blocco_orb = 0;      //vale 1 se ORBGRAND ricostruisce correttamente X
		int undetection_blocco_orb = 0;   //vale 1 se ORBGRAND trova una codeword valida diversa da X
		int abbandono_blocco_orb = 0;     //vale 1 se ORBGRAND raggiunge MAX_QUERY

		//VARIABILI PER SALVARE IL RISULTATO DI CDF-ORBGRAND SUL SINGOLO BLOCCO

		unsigned long long query_blocco_cdf = 0; //numero di query effettuate da CDF-ORBGRAND sul blocco corrente

		int corretto_blocco_cdf = 0;      //vale 1 se CDF-ORBGRAND ricostruisce correttamente X
		int undetection_blocco_cdf = 0;   //vale 1 se CDF-ORBGRAND trova una codeword valida diversa da X
		int abbandono_blocco_cdf = 0;     //vale 1 se CDF-ORBGRAND raggiunge MAX_QUERY


		//CANALE AWGN
		for (i = 0; i < N; i++) {

			//aggiunta del rumore AWGN
			Zgauss = genera_gaussiana_standard();

			//scalatura secondo la \sigma del canale
			RumoreCanale[i] = sigma * Zgauss; //si fa ciò che è stato verificato, ossia Zi con PDF N(0,1) e n_i= \sigma * Z_i, quindi n_i con PDF N(0,\sigma^2)

			//=> per ciascuno dei 128 simboli viene generato un nuovo campione indipendente di rumore 

			//costruisco la Ysoft
			Ysoft[i] = Sbpsk[i] + RumoreCanale[i];

			//dopo aver costruito il singolo valore Ysoft[i], applico la condizione:
			//y_i >= 0 => 0, y_i < 0 => 1
			//=> hard decision => Yhard[i]
			if (Ysoft[i] >= 0) {
				Yhard[i] = 0;
			}
			else {
				Yhard[i] = 1;
			}

			//calcolo degli LLR, viene computato attraverso la relazione:
			// L_i = [ 2 * sqrt(Es) * y_i ] / (\sigma)^2
			LLR[i] = (2.0 * sqrt(Es) * Ysoft[i]) / sigma_quadro;

			//per ricavare l'affidabilità del bit considero il modulo dell'LLR, tramite fabs(), funzione che restituisce il valore assoluto di un double
			Affidabilita[i] = fabs(LLR[i]);

			//nel caso di modulazione standard con Es unitario, la relazione si riduce a LLR = 2 volte Ysoft_i / (\sigma quadro)

			//inizializzazione della permutazione
			pi[i] = i;
		}

		//ordino dal bit meno affidabile al più affidabile tramite la funzione ordina_affidabilita()
		ordina_affidabilita(Affidabilita, pi);

		//CONTROLLO DELLA DECISIONE HARD
		//verifico che ogni elemento di Yhard sia coerente con il segno del corrispondente valore Ysoft

		int controllo_hard = 1;

		for (i = 0; i < N; i++) {

			if (Ysoft[i] >= 0 && Yhard[i] != 0) {

				controllo_hard = 0;
				printf("Errore nella decisione hard alla posizione %d.\n", i);
				break;
			}

			if (Ysoft[i] < 0 && Yhard[i] != 1) {

				controllo_hard = 0;
				printf("Errore nella decisione hard alla posizione %d.\n", i);
				break;
			}
		}

		if (controllo_hard == 1) {

			if (STAMPE_DETTAGLIATE) {
				printf("\nControllo decisione hard superato: Yhard e' coerente con il segno di Ysoft.\n");
			}

			//CONFRONTO TRA CALCOLO DENSO E CALCOLO SPARSO DELLA SINDROME

			calcola_sindrome(Yhard, H, SZ_densa_test);
			calcola_sindrome_sparsa(Yhard, SZ_sparsa_test);

			int sindromi_uguali = 1;

			for (r = 0; r < NH; r++) {

				if (SZ_densa_test[r] != SZ_sparsa_test[r]) {

					sindromi_uguali = 0;
					break;
				}
			}

			if (sindromi_uguali == 1) {

				if (STAMPE_DETTAGLIATE) {
					printf("Controllo sindrome sparsa superato: coincide con il calcolo denso.\n");
				}
			}

			else {

				printf("ERRORE: la sindrome sparsa NON coincide con il calcolo denso.\n");
			}
		}

		//CONTROLLO DEL CALCOLO DELLE AFFIDABILITA'
		//ogni affidabilita' deve coincidere con il valore assoluto del corrispondente LLR

		int controllo_affidabilita = 1;
		double differenza;

		for (i = 0; i < N; i++) {

			differenza = fabs(Affidabilita[i] - fabs(LLR[i]));

			if (differenza > 1e-12) {

				controllo_affidabilita = 0;
				printf("Errore nel calcolo dell'affidabilita' alla posizione %d.\n", i);
				break;
			}
		}

		if (controllo_affidabilita == 1) {

			if (STAMPE_DETTAGLIATE) {
				printf("Controllo affidabilita' superato: Affidabilita[i] = |LLR[i]| per tutti i bit.\n");
			}
		}

		//CONTROLLO DELL'ORDINAMENTO DELLE AFFIDABILITA'
		//dopo ordina_affidabilita(), pi deve indicare i bit in ordine di affidabilita' crescente

		int controllo_ordinamento = 1;

		for (i = 0; i < N - 1; i++) {

			if (Affidabilita[pi[i]] > Affidabilita[pi[i + 1]]) {

				controllo_ordinamento = 0;

				printf("Errore nell'ordinamento tra i ranghi %d e %d.\n",
					i + 1, i + 2);

				break;
			}
		}

		if (controllo_ordinamento == 1) {
			if (STAMPE_DETTAGLIATE) {
				printf("Controllo ordinamento superato: pi ordina i bit dal meno affidabile al piu' affidabile.\n");
			}
		}

		if (STAMPE_DETTAGLIATE) {

			printf("\n===============================================\n");
			printf("TRASMISSIONE BPSK SU CANALE AWGN\n");
			printf("===============================================\n");
		}

		if (STAMPE_DETTAGLIATE) {

			printf("\nX binaria:\n");

			for (i = 0; i < N; i++) {
				printf("%d ", X[i]);

				if ((i + 1) % 16 == 0) {
					printf("\n");
				}
			}
		}

		if (STAMPE_DETTAGLIATE) {

			printf("\nS BPSK:\n");

			for (i = 0; i < N; i++) {
				printf("%+.4f ", Sbpsk[i]);

				if ((i + 1) % 8 == 0) {
					printf("\n");
				}
			}
		}

		if (STAMPE_DETTAGLIATE) {

			printf("\nRumore:\n");

			for (i = 0; i < N; i++) {
				printf("%+.4f ", RumoreCanale[i]);

				if ((i + 1) % 8 == 0) {
					printf("\n");
				}
			}
		}

		if (STAMPE_DETTAGLIATE) {

			printf("\nY soft:\n");

			for (i = 0; i < N; i++) {
				printf("%+.4f ", Ysoft[i]);

				if ((i + 1) % 8 == 0) {
					printf("\n");
				}
			}
		}


		if (STAMPE_DETTAGLIATE) {

			printf("\n===============================================\n");
			printf("HARD DECISION E AFFIDABILITA'\n");
			printf("===============================================\n");
		}

		if (STAMPE_DETTAGLIATE) {

			printf("\nY hard:\n");

			for (i = 0; i < N; i++) {
				printf("%d ", Yhard[i]);

				if ((i + 1) % 16 == 0) {
					printf("\n");
				}
			}
		}

		if (STAMPE_DETTAGLIATE) {

			printf("\nLLR:\n");

			for (i = 0; i < N; i++) {
				printf("%+.4f ", LLR[i]);

				if ((i + 1) % 8 == 0) {
					printf("\n");
				}
			}
		}

		if (STAMPE_DETTAGLIATE) {

			printf("\n|LLR|:\n");

			for (i = 0; i < N; i++) {
				printf("%.4f ", Affidabilita[i]);

				if ((i + 1) % 8 == 0) {
					printf("\n");
				}
			}
		}

		if (STAMPE_DETTAGLIATE) {

			printf("\nPermutazione pi (dal meno affidabile al piu' affidabile):\n");

			for (i = 0; i < N; i++) {
				printf("%d ", pi[i] + 1);

				if ((i + 1) % 16 == 0) {
					printf("\n");
				}
			}

			printf("\n");
		}


		{
			//DECODIFICA ORBGRAND

			//MISURA DEL TEMPO DI DECODIFICA
			clock_t inizio_decoder = clock();

			trovato = 0;

			//NUOVO RISPETTO AL VECCHIO PROGETTO:
			//azzero il contatore prima di iniziare una nuova decodifica
			numero_query = 0;

			//NUOVO ANCORA:
			//all'inizio della decodifica non si è ancora verificato alcun abbandono
			abbandono = 0;

			//ORBGRAND esplora automaticamente i Logistic Weight in ordine crescente
			for (wL = 0; wL <= wL_max; wL++) {

				if (STAMPE_DETTAGLIATE) {

					printf("\n=====================================\n");
					printf("LOGISTIC WEIGHT wL = %d\n", wL);
					printf("=====================================\n");
				}

				//ogni nuovo Logistic Weight deve partire da un pattern di rango nullo
				azzera_vettore(Erank, N);

				//GENERAZIONE E PROVA DEI PATTERN LOGISTICI CON PESO wL

				//MODIFICA RISPETTO AL VECCHIO PROGETTO:
				//Nel progetto toy veniva passata la matrice L:
				//genera_pattern_logistici(wL, 1, Erank, Ecan, pi, Yhard, L, Z, SZ, Xcan, Efin, &trovato);
				//Nel nuovo progetto LDPC viene invece passata la matrice originale H.
				//L viene utilizzata come matrice di lavoro ridotta per la costruzione di G,
				//Mentre H rappresenta la matrice originale di controllo di parità che definisce il codice.
				//Di conseguenza, le candidate costruite da ORBGRAND vengono verificate tramite Z * H^T = 0.

				genera_pattern_logistici(wL, 1, Erank, Ecan, pi, Yhard, H, Z, SZ, Xcan, Efin, &trovato);

				//questa chiamata, ad ogni iterazione, presenta:
				//peso_rimanente = wL e rango_iniziale = 1
				//=> per wL=0 viene provato il pattern nullo;
				//se fallisce, il ciclo passa a wL=1, poi 2 eccetera

				//appena viene trovata una codeword valida il decoder si arresta
				//if (trovato == 1) {
					//break;
				//}
					
				//il decoder termina se ha trovato una codeword valida
				//oppure se ha raggiunto il numero massimo di query
				if (trovato == 1 || abbandono == 1) {
					break;
				}
			}

			clock_t fine_decoder = clock();

			double tempo_decoder =
				(double)(fine_decoder - inizio_decoder) / CLOCKS_PER_SEC;

			//BLOCCO FINALE DI STAMPA E VERIFICA DEL DECODER

			printf("\n=====================================\n");
			printf("RISULTATO FINALE DEL DECODER\n");
			printf("=====================================\n");

			printf("Tempo di decodifica: %.6f secondi\n", tempo_decoder);

			if (trovato == 1) {

				if (numero_query > query_massime_successo_orb) {
					query_massime_successo_orb = numero_query;
				}

				//MODIFICA RISPETTO AL VECCHIO PROGETTO:
				//nel progetto toy veniva utilizzata stampa_vettore().
				//Nel nuovo progetto, con N=128, i vettori vengono stampati direttamente
				//disponendo gli elementi su più righe per rendere l'output più leggibile.

				if (STAMPE_DETTAGLIATE) {

					printf("\nCodeword stimata Xcan:\n");
					for (i = 0; i < N; i++) {
						printf("%d ", Xcan[i]);

						if ((i + 1) % 16 == 0) {
							printf("\n");
						}
					}

					printf("\nPattern trovato Ecan finale (Efin):\n");
					for (i = 0; i < N; i++) {
						printf("%d ", Efin[i]);

						if ((i + 1) % 16 == 0) {
							printf("\n");
						}
					}

					printf("\nSindrome finale SZ:\n");
					for (i = 0; i < NH; i++) {
						printf("%d ", SZ[i]);

						if ((i + 1) % 16 == 0) {
							printf("\n");
						}
					}
				}

				//MODIFICA RISPETTO AL VECCHIO PROGETTO:
				//la dicitura GRAND viene aggiornata in ORBGRAND,
				//poiché in questo progetto la ricerca è ordinata tramite le affidabilità.

				printf("\nORBGRAND ha trovato una parola valida del codice.\n");

				printf("Numero totale di query effettuate: %llu\n", numero_query);
				printf("Logistic Weight finale raggiunto: %d\n", wL);

				int uguali = 1;

				for (i = 0; i < N; i++) {

					if (X[i] != Xcan[i]) {
						uguali = 0;
						break;
					}
				}

				if (uguali == 1) {

					numero_corretti_orb++;
					corretto_blocco_orb = 1;

					printf("Decodifica corretta: Xcan coincide con la sequenza trasmessa X.\n");
				}
				else {

					numero_undetection_orb++;
					undetection_blocco_orb = 1;

					printf("ORBGRAND ha trovato una parola valida del codice, ma diversa dalla sequenza trasmessa X.\n");
				}
			}
			else {

				if (abbandono == 1) {

					numero_abbandoni_orb++;
					abbandono_blocco_orb = 1;

					printf("\nORBGRAND ha raggiunto il numero massimo di query.\n");
					printf("Numero totale di query effettuate: %llu\n", numero_query);
					printf("Limite massimo di query: %llu\n", (unsigned long long)MAX_QUERY);
					printf("Decodifica interrotta: il blocco viene classificato come ABBANDONO.\n");
				}

				else {

					printf("ORBGRAND non ha trovato una parola valida tra i pattern testati.\n");
				}
			}

			query_blocco_orb = numero_query;
			query_totali_orb += numero_query;
		}
		
		
		//DECODIFICA CDF-ORBGRAND
		//Il pattern nullo viene provato separatamente; se non produce una codeword valida, la ricerca continua estraendo dalla min-heap il pattern con peso CDF minimo

		//MISURA DEL TEMPO DI DECODIFICA
		clock_t inizio_decoder = clock();

		NodoHeapCDF nodo_iniziale_decoder_cdf = { 0 };
		NodoHeapCDF nodo_corrente_decoder_cdf;

		//all'inizio di ogni nuova decodifica la min-heap deve essere vuota
		heap_decoder_cdf.dimensione = 0;

		trovato = 0;

		//azzero il contatore prima di iniziare una nuova decodifica
		numero_query = 0;

		//all'inizio della decodifica non si è ancora verificato alcun abbandono
		abbandono = 0;


		//PRIMA QUERY: PATTERN NULLO
		//Il pattern nullo ha peso CDF nullo e deve sempre essere il primo tentativo del decoder.

		azzera_vettore(Ecan, N);

		prova_pattern(Yhard, H, Ecan, Z, SZ, Xcan, Efin, &trovato);


		//Se il pattern nullo non produce una codeword valida, la min-heap viene inizializzata con il primo pattern non nullo, corrispondente al solo rango di affidabilità 1.
		if (trovato == 0 && abbandono == 0) {

			inserisci_rango_cdf(&nodo_iniziale_decoder_cdf, 1);

			nodo_iniziale_decoder_cdf.peso_cdf = pesi_cdf[0];
			nodo_iniziale_decoder_cdf.rango_massimo = 1;

			inserisci_heap(&heap_decoder_cdf, nodo_iniziale_decoder_cdf);
		}


		//RICERCA CDF
		//Ad ogni iterazione viene estratto il pattern avente peso CDF minimo tra quelli presenti nella heap.
		//Il pattern viene trasformato dalle posizioni di rango alle posizioni fisiche tramite pi e successivamente verificato mediante sindrome.

		while (trovato == 0 &&
			abbandono == 0 &&
			estrai_minimo_heap(&heap_decoder_cdf, &nodo_corrente_decoder_cdf) == 1) {

			//Il pattern espresso nel dominio dei ranghi viene riportato nelle effettive posizioni della sequenza ricevuta.
			mappa_pattern_cdf(nodo_corrente_decoder_cdf, Ecan, pi);

			//La parola candidata viene costruita e verificata mediante sindrome.
			prova_pattern(Yhard, H, Ecan, Z, SZ, Xcan, Efin, &trovato);

			//Se dopo questa query è stato raggiunto MAX_QUERY senza trovare una codeword valida, il blocco viene classificato come abbandono.
			if (trovato == 0 && numero_query >= MAX_QUERY) {

				abbandono = 1;
				break;
			}

			//Se il pattern corrente non produce una codeword valida, vengono generati i nuovi pattern candidati da inserire nella min-heap.
			if (trovato == 0) {

				genera_nuovi_pattern_cdf(nodo_corrente_decoder_cdf, pesi_cdf, &heap_decoder_cdf);
			}
		}

		clock_t fine_decoder = clock();

		double tempo_decoder =
			(double)(fine_decoder - inizio_decoder) / CLOCKS_PER_SEC;

		//BLOCCO FINALE DI STAMPA E VERIFICA DEL DECODER

		printf("\n=====================================\n");
		printf("RISULTATO FINALE DEL DECODER\n");
		printf("=====================================\n");

		printf("Tempo di decodifica: %.6f secondi\n", tempo_decoder);

		if (trovato == 1) {

			if (numero_query > query_massime_successo) {
				query_massime_successo = numero_query;
			}

			//MODIFICA RISPETTO AL VECCHIO PROGETTO:
			//nel progetto toy veniva utilizzata stampa_vettore().
			//Nel nuovo progetto, con N=128, i vettori vengono stampati direttamente
			//disponendo gli elementi su più righe per rendere l'output più leggibile.

			if (STAMPE_DETTAGLIATE) {

				printf("\nCodeword stimata Xcan:\n");
				for (i = 0; i < N; i++) {
					printf("%d ", Xcan[i]);

					if ((i + 1) % 16 == 0) {
						printf("\n");
					}
				}

				printf("\nPattern trovato Ecan finale (Efin):\n");
				for (i = 0; i < N; i++) {
					printf("%d ", Efin[i]);

					if ((i + 1) % 16 == 0) {
						printf("\n");
					}
				}

				printf("\nSindrome finale SZ:\n");
				for (i = 0; i < NH; i++) {
					printf("%d ", SZ[i]);

					if ((i + 1) % 16 == 0) {
						printf("\n");
					}
				}
			}

			printf("\nCDF-ORBGRAND ha trovato una parola valida del codice.\n");

			printf("Numero totale di query effettuate: %llu\n", numero_query);

			int uguali = 1;

			for (i = 0; i < N; i++) {

				if (X[i] != Xcan[i]) {
					uguali = 0;
					break;
				}
			}

			if (uguali == 1) {

				numero_corretti++;
				corretto_blocco_cdf = 1;

				printf("Decodifica corretta: Xcan coincide con la sequenza trasmessa X.\n");
			}
			else {

				numero_undetection++;
				undetection_blocco_cdf = 1;

				printf("CDF-ORBGRAND ha trovato una parola valida del codice, ma diversa dalla sequenza trasmessa X.\n");
			}
		}
		else {

			if (abbandono == 1) {

				numero_abbandoni++;
				abbandono_blocco_cdf = 1;

				printf("\nCDF-ORBGRAND ha raggiunto il numero massimo di query.\n");
				printf("Numero totale di query effettuate: %llu\n", numero_query);
				printf("Limite massimo di query: %llu\n", (unsigned long long)MAX_QUERY);
				printf("Decodifica interrotta: il blocco viene classificato come ABBANDONO.\n");
			}

			else {

				printf("CDF-ORBGRAND non ha trovato una parola valida tra i pattern testati.\n");
			}
		}

		query_blocco_cdf = numero_query;
		query_totali += numero_query;

		//CONFRONTO BLOCCO PER BLOCCO TRA ORBGRAND E CDF-ORBGRAND
		//Entrambi i decoder hanno lavorato sullo stesso messaggio trasmesso e sulla stessa realizzazione di rumore => il confronto delle query è effettuato sullo stesso identico blocco

		printf("\n=====================================\n");
		printf("CONFRONTO SUL BLOCCO %d\n", numero_blocchi);
		printf("=====================================\n");

		printf("ORBGRAND: ");

		if (corretto_blocco_orb == 1) {
			printf("CORRETTO");
		}
		else if (undetection_blocco_orb == 1) {
			printf("UNDETECTION");
		}
		else if (abbandono_blocco_orb == 1) {
			printf("ABBANDONO");
		}

		printf(" - Query: %llu\n", query_blocco_orb);


		printf("CDF-ORBGRAND: ");

		if (corretto_blocco_cdf == 1) {
			printf("CORRETTO");
		}
		else if (undetection_blocco_cdf == 1) {
			printf("UNDETECTION");
		}
		else if (abbandono_blocco_cdf == 1) {
			printf("ABBANDONO");
		}

		printf(" - Query: %llu\n", query_blocco_cdf);


		//CASO DI PARTICOLARE INTERESSE PER IL CONFRONTO
		//ORBGRAND raggiunge MAX_QUERY mentre CDF-ORBGRAND riesce a recuperare correttamente lo stesso blocco
		if (abbandono_blocco_orb == 1 && corretto_blocco_cdf == 1) {

			recuperi_cdf++;
			
			printf("=> CDF-ORBGRAND RECUPERA UN BLOCCO ABBANDONATO DA ORBGRAND.\n");
		}

		//VICEVERSA
		if (abbandono_blocco_cdf == 1 && corretto_blocco_orb == 1) {

			recuperi_orb++;

			printf("=> ORBGRAND RECUPERA UN BLOCCO ABBANDONATO DA CDF-ORBGRAND.\n");
		}

		//CONFRONTO DEL NUMERO DI QUERY QUANDO ENTRAMBI I DECODER CORREGGONO LO STESSO BLOCCO
		if (corretto_blocco_orb == 1 && corretto_blocco_cdf == 1) {

			if (query_blocco_cdf < query_blocco_orb) {

				cdf_meno_query++;
			}
			else if (query_blocco_orb < query_blocco_cdf) {

				orb_meno_query++;
			}
			else {

				query_uguali++;
			}
		}


	} //fine ciclo sui blocchi

	//RIEPILOGO DEL TEST SU PIU' BLOCCHI

	int numero_errori = numero_undetection + numero_abbandoni; //CDF
	int numero_errori_orb = numero_undetection_orb + numero_abbandoni_orb; //ORB

	//CALCOLO DELLA BLOCK ERROR RATE
	//numero_errori rappresenta il numero complessivo di blocchi errati osservati, ottenuto come somma tra undetection e abbandoni
	//al termine del while, in condizioni normali, numero_errori = TARGET_ERRORI = 100
	//numero_blocchi rappresenta, invece, il numero totale M di blocchi che è stato necessario simulare per osservare tali 100 errori
	//=> BLER = numero_errori / numero_blocchi = 100 / M

	double BLER = (double)numero_errori / numero_blocchi; //CDF
	double BLER_orb = (double)numero_errori_orb / numero_blocchi; //ORB

	double query_medie =
		(double)query_totali / numero_blocchi; //CDF

	double query_medie_orb =
		(double)query_totali_orb / numero_blocchi; //ORB

	printf("\n\n=========================================\n");
	printf("RIEPILOGO SIMULAZIONE\n");
	printf("=========================================\n");

	printf("Blocchi totali comuni: %d\n", numero_blocchi);

	printf("\nORBGRAND\n");
	printf("Decodifiche corrette: %d\n", numero_corretti_orb);
	printf("Undetection: %d\n", numero_undetection_orb);
	printf("Abbandoni: %d\n", numero_abbandoni_orb);
	printf("Errori totali: %d\n", numero_errori_orb);
	printf("BLER: %.8f\n", BLER_orb);
	printf("Query totali: %llu\n", query_totali_orb);
	printf("Query medie per blocco: %.2f\n", query_medie_orb);
	printf("Massimo numero di query tra i successi: %llu\n", query_massime_successo_orb);

	printf("\nCDF-ORBGRAND\n");
	printf("Decodifiche corrette: %d\n", numero_corretti);
	printf("Undetection: %d\n", numero_undetection);
	printf("Abbandoni: %d\n", numero_abbandoni);
	printf("Errori totali: %d\n", numero_errori);
	printf("BLER: %.8f\n", BLER);
	printf("Query totali: %llu\n", query_totali);
	printf("Query medie per blocco: %.2f\n", query_medie);
	printf("Massimo numero di query tra i successi: %llu\n", query_massime_successo);

	printf("\nCONFRONTO DIRETTO\n");
	printf("Blocchi recuperati da CDF-ORBGRAND rispetto a ORBGRAND: %d\n", recuperi_cdf);
	printf("Blocchi recuperati da ORBGRAND rispetto a CDF-ORBGRAND: %d\n", recuperi_orb);

	printf("Blocchi corretti da entrambi con meno query per CDF-ORBGRAND: %d\n", cdf_meno_query);
	printf("Blocchi corretti da entrambi con meno query per ORBGRAND: %d\n", orb_meno_query);
	printf("Blocchi corretti da entrambi con lo stesso numero di query: %d\n", query_uguali);

	free(heap_decoder_cdf.nodi);
	free(M);
	free(G);

	return 0;
}