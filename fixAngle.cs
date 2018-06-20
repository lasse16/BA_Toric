using System;
using UnityEngine;
/**
 * A representation of an angle, always under 360 degrees and with a easy toRad function.
 * 
 * @author Lasse Haffke 
 * @version 20.06.18
 */

public class FixAngle
{
	private float _angle;
	
	public FixAngle(float angle)
	{
		_angle = angle % 360;
	}
	
	public float angle()
	{
		return _angle;
	}
	
	public float toRad(){
		return _angle * Mathf.Deg2Rad;
	}
}

