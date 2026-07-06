#pragma once

#include <AudioToolbox/AudioToolbox.h>
#include "SoundChainPlatform.hpp"

class CoreAudioSCP : public SoundChainPlatform {
public:
	CoreAudioSCP() {}
	~CoreAudioSCP() {}

	static OSStatus render_callback(void *refCon, 
			AudioUnitRenderActionFlags *flags, 
			const AudioTimeStamp *ts, 
			UInt32 bus, 
			UInt32 frames, 
			AudioBufferList *ioData) {

		CoreAudioSCP *scPlatform = (CoreAudioSCP *)refCon;
		float *out = (float *)ioData->mBuffers[0].mData;

		scPlatform->FillBuffer(out, (int)frames);

		return noErr;
	};

private:
	AudioUnit _unit = {};

	void Setup() override {
		AudioComponentDescription desc = {0};
		desc.componentType = kAudioUnitType_Output;
		desc.componentSubType = kAudioUnitSubType_DefaultOutput;
		desc.componentManufacturer = kAudioUnitManufacturer_Apple;

		AudioComponent comp = AudioComponentFindNext(nullptr, &desc);
		if (!comp) exit(EXIT_FAILURE);

		OSStatus err = AudioComponentInstanceNew(comp, &_unit);
		if (err != noErr) exit(EXIT_FAILURE);

		AURenderCallbackStruct cb = {0};
		cb.inputProc = render_callback;
		cb.inputProcRefCon = this;
		
		err = AudioUnitSetProperty(
			_unit,
			kAudioUnitProperty_SetRenderCallback,
			kAudioUnitScope_Input,
			0,
			&cb,
			sizeof(cb)
		);
		if (err != noErr) exit(EXIT_FAILURE);

		AudioStreamBasicDescription streamDesc = {};
		streamDesc.mSampleRate = GetSoundChainSettings().SampleRate;
		streamDesc.mFormatID = kAudioFormatLinearPCM;
		streamDesc.mFormatFlags = kAudioFormatFlagIsFloat;
		streamDesc.mBitsPerChannel = 32;
		streamDesc.mChannelsPerFrame = (UInt32)GetSoundChainSettings().Channels;
		streamDesc.mFramesPerPacket = 1;
		streamDesc.mBytesPerFrame = 4 * streamDesc.mChannelsPerFrame;
		streamDesc.mBytesPerPacket = streamDesc.mBytesPerFrame * streamDesc.mFramesPerPacket;

		// Non-interleaved: set on Input scope (commonly for output units)
		err = AudioUnitSetProperty(
			_unit,
			kAudioUnitProperty_StreamFormat,
			kAudioUnitScope_Input,
			0,
			&streamDesc,
			sizeof(streamDesc)
		);
		if (err != noErr) exit(EXIT_FAILURE);

		err = AudioUnitInitialize(_unit);
		if (err != noErr) exit(EXIT_FAILURE);
	}

	void Start() override {
		AudioOutputUnitStart(_unit);
	}

	void End() override {
		AudioOutputUnitStop(_unit);
		AudioUnitUninitialize(_unit);
		AudioComponentInstanceDispose(_unit);
	}
};