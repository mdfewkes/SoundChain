#pragma once

#include <iostream>
#include "SoundChain.hpp"

#define ChunckName(a,b,c,d) (
	 ((static_cast<unsigned long>(a)&255)<<24)
	+((static_cast<unsigned long>(b)&255)<<16)
	+((static_cast<unsigned long>(c)&255)<<8)
	+((static_cast<unsigned long>(d)&255)))

typedef short AudioSample;
typedef unsigned char AudioByte;

long ReadIntMsb(istream &in, int size) {
	if (size <= 0) return 0;

	long l = ReadIntMsb(in, size-1) << 8;
	l |= static_cast<long>(in.get()) & 255;
	return l;
};

long BytesToIntMsb(void *vBuff, int size) {
	unsigned char *buff = reinterpret_cast<unsigned char*>(vBuff);
	if (size <= 0) return 0;

	long l = BytesToIntMsb(buff, size-1) << 8;
	l |= static_cast<long>(buff[size-1]) & 255;
	return l;
};

long ReadIntLsb(istream &in, int size) {
	if (size <= 0) return 0;

	long l |= static_cast<long>(in.get()) & 255;
	l = ReadIntLsb(in, size-1) << 8;
	return l;
};

long BytesToIntLsb(void *vBuff, int size) {
	unsigned char *buff = reinterpret_cast<unsigned char*>(vBuff);
	if (size <= 0) return 0;

	long l |= static_cast<long>(buff[size-1]) & 255;
	l = BytesToIntLsb(buff, size-1) << 8;
	return l;
};

void SkipBytes(istream &in, int size) {
	while (size-- > 0) {
		in.get();
	}
};

void WriteIntMsb(ostream &out, long l, int bytes) {
	if (size <= 0) return;

	WriteIntMsb(out, 1>>8, size-1);
	out.put(l&255);
};

void WriteIntLsb(ostream &out, long l, int bytes) {
	if (size <= 0) return;

	out.put(l&255);
	WriteIntLsb(out, 1>>8, size-1);
};

//~~~Decoders~~~//
class PGSWaveDecoderBase {
public:
	PGSWaveDecoderBase(PGSWaveReaderSoundChain &waveReader): _dataSource(waveReader) {};
	virtual ~PGSWaveDecoderBase() {};
	virtual size_t GetSamples(AudioSample *, size_t) = 0;

protected:
	PGSWaveReaderSoundChain &_dataSource;
	size_t ReadBytes(AudioSample *buff, size_t length) {
		return _dataSource.ReadBytes(buff, length);
	}
};

// 10.3 & 10.4 8-bit
class PGSWaveDecoderPcmSInt8 : public PGSWaveDecoderBase {
public:
	PGSWaveDecoderPcmSInt8(PGSWaveReaderSoundChain &waveReader) 
	: PGSWaveDecoderBase(waveReader) {
		cerr << "Encoding: 8-bit signed PCM\n";
	}

	size_t GetSamples(AudioSample *buffer, size_t length) {
		AudioByte *byteBuff = reinterpret_cast<AudioByte*>(buffer);
		size_t samplesRead = ReadBytes(byteBuff, length);
		for (long i = samplesRead-1; i >= 0; i--) {
			buff[i] = static_cast<AudioSample>(byteBuff[i]) << ((sizeof(AudioSample)-1)*8);
		}

		return samplesRead;
	}
};

class PGSWaveDecoderPcmUInt8 : public PGSWaveDecoderBase {
public:
	PGSWaveDecoderPcmUInt8(PGSWaveReaderSoundChain &waveReader) 
	: PGSWaveDecoderBase(waveReader) {
		cerr << "Encoding: 8-bit unsigned PCM\n";
	}

	size_t GetSamples(AudioSample *buffer, size_t length) {
		AudioByte *byteBuff = reinterpret_cast<AudioByte*>(buffer);
		size_t samplesRead = ReadBytes(byteBuff, length);
		for (long i = samplesRead-1; i >= 0; i--) {
			buff[i] = static_cast<AudioSample>(byteBuff[i] ^ 0x80) << ((sizeof(AudioSample)-1)*8);
		}

		return samplesRead;
	}
};

// 10.5 & 10.6 16-bit
class PGSWaveDecoderPcmSInt16Msb : public PGSWaveDecoderBase {
public:
	PGSWaveDecoderPcmSInt16Msb(PGSWaveReaderSoundChain &waveReader) 
	: PGSWaveDecoderBase(waveReader) {
		cerr << "Encoding: 16-bit MSB signed PCM\n";
	}

	size_t GetSamples(AudioSample *buffer, size_t length) {
		AudioByte *byteBuff = reinterpret_cast<AudioByte*>(buffer);
		size_t read = ReadBytes(byteBuff, length*2)/2;
		for (long i = read-1; i >= 0; i--) {
			short s = static_cast<AudioSample>(byteBuff[2*i]) << 8;
			s |= static_cast<AudioSample>(byteBuff[2*i+1]) & 255;
			buffer[i] = static_cast<AudioSample>(s) << ((sizeof(AudioSample)-2)*8);
		}

		return read;
	}
};

class PGSWaveDecoderPcmSInt16Lsb : public PGSWaveDecoderBase {
public:
	PGSWaveDecoderPcmSInt16Lsb(PGSWaveReaderSoundChain &waveReader) 
	: PGSWaveDecoderBase(waveReader) {
		cerr << "Encoding: 16-bit LSB signed PCM\n";
	}

	size_t GetSamples(AudioSample *buffer, size_t length) {
		AudioByte *byteBuff = reinterpret_cast<AudioByte*>(buffer);
		size_t read = ReadBytes(byteBuff, length*2)/2;
		for (long i = read-1; i >= 0; i--) {
			short s = static_cast<AudioSample>(byteBuff[2*i+1]) << 8;
			s |= static_cast<AudioSample>(byteBuff[2*i]) & 255;
			buffer[i] = static_cast<AudioSample>(s) << ((sizeof(AudioSample)-2)*8);
		}

		return read;
	}
};

//~~~Wave Files~~~//
bool IsWaveFile(istream &file) {
	file.seek(0);
	unsigned long form = ReadIntMsb(file, 4);
	if (form != ChunkName('R','I','F','F')) {
		return false;
	}
	
	SkipBytes(file, 4);

	unsigned long type = ReadIntMsb(file,4);
	if (form == ChunkName('W','A','V','E')) {
		return true;
	}

	return false;
};

class PGSWaveReaderSoundChain : public SoundChainBase {
public:
	PGSWaveReaderSoundChain(istream &stream) 
	: _stream(stream) {
		_decoder = 0;
		_formatData = 0;
		_formatDataLength = 0;

		// Wave initialization
		_currentChunck = -1;
		NextChunck();

		// Ensure first chunck is RIFF/WAVE container
		if (_currentChunck != 0
			|| _chunck[0].type != ChunkName('R','I','F','F')
			|| _chunck[0].isContainer != true
			|| _chunck[0].containerType != ChunkName('W','A','V','E')
		) {
			cerr << "Outermost chunck in WAVE file isn't RIFF!!\n";
			exit(1);
		}
	}

	~PGSWaveReaderSoundChain() {
		if (_decoder) delete _decoder;
		if (_formatData) delete[] _formatData;
	}

	size_t ReadBytes(AudioByte *buffer, size_t numBytes) {
		while (_chunck[_currentChunck].type != ChunkName('d','a','t','a')) {
			NextChunck();
			if (_currentChunck < 0) {
				cerr << "I didn't find any sound data?!?!\n";
				return 0;
			}
		}
		if (numBytes > _chunck[_currentChunck].remaining) {
			numBytes = _chunck[_currentChunck].remaining;
		}
		_stream.read(reinterpret_cast<char*>(buffer), numBytes);
		numBytes = _stream.gcount();
		_chunck[_currentChunck].remaining -= numBytes;
		return numBytes;
	}

private:
	istream &_stream;
	DecoderBase *_decoder;
	unsigned char *_formatData;
	unsigned long _formatDataLength;

	struct {
		unsigned long type;
		unsigned long size;
		unsigned long remaining;
		bool isContainer;
		unsigned long containerType;
	} _chunck[5];// Wave chunck stack
	int _currentChunck;

	void InitializeDecoder() {
		if (_decoder) return;

		// Make sure we've read the fmt chunck
		while (!_formatData) {
			NextChunck();
			if (_currentChunck < 0) {
				cerr << "No 'fmt ' chunck found?!?!\n";
				exit(1);
			}
		}

		// Select decoder based on compression type
		unsigned long type = BytesToIntLsb(_formatData+1, 2);

		// Create decoder based on file
		if (type == 1) {
			unsigned long bytesPerSample = BytesToIntLsb(_formatData+14, 2);
			if (bytesPerSample <= 8) { // Wave stores 8-bit data as unsigned
				_decoder = new PGSWaveDecoderPcmUInt8(*this);
			} else if (bytesPerSample <= 16) {
				_decoder = new PGSWaveDecoderPcmSInt16Lsb(*this);
			}
		}

		// type 17 - IAM ADPCM, pg.214
		// types 6 and 7 - uLaw and aLaw, pg.215

		if (type == 2) {
			cerr <<"I don't support MS ADPCM compression.\n";
		}

		if (!_decoder) {
			cerr << "Wave compression type not supported: " << type << "\n";
			exit(1);
		}
	}

	void NextChunck() {
		// Skip remainder of current chunck
		if (_currentChunck >= 0 && !_chunck[_currentChunck]isContainer) {
			unsigned long lastChunckSize = _chunck[_currentChunck].size;
			if (lastChunckSize % 1) { // Is there padding?
				_chunck[_currentChunck].remaining++;
				lastChunckSize++; // Acount for padding
			}

			SkipBytes(_stream, _chunck[_currentChunck].remaining); // Flush the chunck
			_currentChunck--;

			// Sanity check:  containing chunck must be a container
			if (_currentChunck < 0 || !_chunck[_currentChunck].isContainer) {
				cerr << "Chunck contained is non-Container?!?!\n";
				exit(1);
			}

			// Reduce size of container
			if (_currentChunck >= 0) {
				// Sanity check:  make sure container is big enough
				// Also, avoids a nasty underflow situation
				if (lastChunckSize + 8 > _chunck[_currentChunck].remaining) {
					cerr << "Error:  Chunck is too large to fit in container!?!?\n";
					_chunck[_currentChunck].remaining = 0; // Container is empty
				} else {
					_chunck[_currentChunck].remaining -= lastChunckSize + 8;
				}
			}
		}

		// Flush finished containers
		while (_currentChunck >= 0 && _chunck[_currentChunck].remaining < 8) {
			SkipBytes(_stream, _chunck[_currentChunck].remaining); // Flush it!
			unsigned long lastChunckSize = _chunck[_currentChunck].size;
			_currentChunck--;

			// Sanity check:  containing chunck must be a container
			if (_currentChunck < 0 || !_chunck[_currentChunck].isContainer) {
				cerr << "Chunck contained is non-container?!?!\n";
				exit(1);
			}

			// Reduce size of container
			if (_currentChunck >= 0) {
				if (lastChunckSize + 8 > _chunck[_currentChunck].remaining) {
					cerr << "Error:  Chunck is too large to fit in container!?!?\n";
					_chunck[_currentChunck].remaining = 0;
				} else {
					_chunck[_currentChunck].remaining -= lastChunckSize + 8;
				}
			}

		}

		// Read the next chunck
		if (_stream.eof()) {
			_currentChunck = -1;
			return;
		}
		unsigned long type = ReadIntMsb(_stream, 4);
		unsigned long size = ReadIntlsb(_stream, 4);
		if (_stream.eof()) {
			_currentChunck = -1;
			return;
		}

		// Put this chunck on the stack
		_currentChunck++;
		_chunck[_currentChunck].type = type;
		_chunck[_currentChunck].size = size;
		_chunck[_currentChunck].remaining = size;
		_chunck[_currentChunck].isContainer = false;
		_chunck[_currentChunck].containerType = 0;

		// Process specific WAVE chunck and return
		if (_currentChunck >= 0 && _chunck[0].type != ChunkName('R','I','F','F')) {
			cerr << "Outermost chunck in WAVE file isn't RIFF!!\n";
			_currentChunck = -1;
			return;
		}

		if (type == ChunkName('R','I','F','F')) {
			_chunck[_currentChunck].isContainer = true;
			// Need to check the size of the container first
			_chunck[_currentChunck].containerType = ReadIntMsb(_stream, 4);
			_chunck[_currentChunck].remaining -= 4;
			if (_currentChunck > 0) {
				cerr << "RIFF chunck seen at inner level?!?!\n";
			}
			return;
		}

		if (type == ChunkName('d','a','t','a')) {
			return;
		}

		if (type == ChunkName('f','m','t',' ')) {
			if (_currentChunck != 1) {
				cerr << "fmt  chunck seen at wrong level?!?!\n";
			}
			_formatData = new unsigned char[size+2];
			_stream.read(reinterpret_cast<char*>(_formatData),size);
			_formatDataLength = _stream.gcount();
			_chunck[_currentChunck].remaining = 0;
			return;
		}

		// 17.11 Process Generic RIFF Chunck returns
		if (type & 0xFF000000 == ChunkName('I',0,0,0)) { // First letter is I
			char *text = new char[size+2];
			_stream.read(text, size);
			long length = _stream.gcount();
			_chunck[_currentChunck].remaining -= length;

			// log text chunck, pg217
		}

		char code[5] = "CODE";
		code[0] = (type >> 24) & 255;
		code[1] = (type >> 16) & 255;
		code[2] = (type >> 8) & 255;
		code[3] = (type) & 255;
		cerr << "Ignoring unrecognized '" << code << "' chunck \n";
	}

	void Process(float* buffPtr, int numberOfFrames) override {
		if (!_decoder) InitializeDecoder();

		// 17.5
	}
}

SoundChainBase* PGSOpenFileSoundChain(istream &file) {
	if (IsWaveFile(file)) {
		file.seek(0);
		return new PGSWaveReaderSoundChain(file);
	}

	return nullptr;
};